# 03_rf_train.R - Random forest: split, tuning (out-of-bag), benchmark against
# the published model and the empirical models
source("R/00_setup.R")
library(rsample)

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
df  <- dat$model_df
pv  <- dat$predictors
fml <- as.formula(paste(TARGET, "~", paste(pv, collapse = " + ")))

# ---- 70/30 split stratified on the response ----------------------------------
set.seed(SEED)
sp     <- initial_split(df, prop = 0.70, strata = TARGET)
train  <- training(sp)
test   <- testing(sp)
cat(sprintf("Split: train=%d  test=%d\n", nrow(train), nrow(test)))

n_threads <- max(1L, parallel::detectCores() - 1L)

# ---- Tuning with ranger (out-of-bag RMSE) -------------------------------------
grid <- expand.grid(
  mtry            = c(3, 4, 5, 6, 8, 10),
  min.node.size   = c(5, 10, 17, 25, 50),
  replace         = c(TRUE, FALSE),
  sample.fraction = c(0.5, 0.632, 0.8, 1.0),
  stringsAsFactors = FALSE
)
# sample.fraction = 1 with replace = FALSE is a degenerate bootstrap: drop it
grid <- grid[!(grid$replace == FALSE & grid$sample.fraction == 1.0), ]
cat("Tuning configurations:", nrow(grid), "\n")

oob <- numeric(nrow(grid))
for (i in seq_len(nrow(grid))) {
  m <- ranger(fml, data = train, num.trees = 500,
              mtry            = grid$mtry[i],
              min.node.size   = grid$min.node.size[i],
              replace         = grid$replace[i],
              sample.fraction = grid$sample.fraction[i],
              num.threads = n_threads, seed = SEED)
  oob[i] <- sqrt(m$prediction.error)        # out-of-bag RMSE
}
grid$oob_rmse <- oob
grid <- grid[order(grid$oob_rmse), ]
best <- grid[1, ]
cat("\nBest configuration (out-of-bag):\n"); print(best, row.names = FALSE)
write.csv(grid, file.path(PATHS$tables, "T3_rf_tuning.csv"), row.names = FALSE)

# ---- num.trees from the out-of-bag curve -------------------------------------
ntrees <- seq(100, 1500, by = 100)
oob_nt <- sapply(ntrees, function(nt) {
  m <- ranger(fml, data = train, num.trees = nt,
              mtry = best$mtry, min.node.size = best$min.node.size,
              replace = best$replace, sample.fraction = best$sample.fraction,
              num.threads = n_threads, seed = SEED)
  sqrt(m$prediction.error)
})
best_ntree <- ntrees[which.min(oob_nt)]
cat("Best num.trees (out-of-bag):", best_ntree, "\n")
write.csv(data.frame(ntree = ntrees, oob_rmse = oob_nt),
          file.path(PATHS$tables, "T3b_rf_ntrees.csv"), row.names = FALSE)

# ---- Final model: ranger (predictive) ----------------------------------------
rf_ranger <- ranger(fml, data = train, num.trees = best_ntree,
                     mtry = best$mtry, min.node.size = best$min.node.size,
                     replace = best$replace, sample.fraction = best$sample.fraction,
                     importance = "permutation", num.threads = n_threads, seed = SEED)

# ---- Final model: randomForest (e2tree back end) -----------------------------
# nodesize ~ min.node.size; sampling consistent with replace/sample.fraction
ss <- if (best$replace) nrow(train) else floor(best$sample.fraction * nrow(train))
set.seed(SEED)
rf <- randomForest(fml, data = train, ntree = best_ntree,
                    mtry = best$mtry, nodesize = best$min.node.size,
                    replace = best$replace, sampsize = ss,
                    importance = TRUE)

# ---- Evaluation on the test set ----------------------------------------------
pred_ranger <- predict(rf_ranger, test)$predictions
pred_rf     <- predict(rf, test)

res <- rbind(
  metrics_tbl(test[[TARGET]], pred_ranger, "RF (ranger, final predictive model)"),
  metrics_tbl(test[[TARGET]], pred_rf,     "RF (randomForest, e2tree back end)")
)

# Comparison with the empirical models (predictions already in a column, in cm)
for (em in EMPIRICAL_COLS) if (em %in% names(test))
  res <- rbind(res, metrics_tbl(test[[TARGET]], test[[em]], em))

# Published reference (Ebrahimi et al. 2024, Table 4, test set, different n)
res_ref <- data.frame(model = "RF published (Table 4, test set, n=1843)",
                       RMSE = 82.46, MAE = 55.20, R2 = NA, Bias = NA)
res <- rbind(res, res_ref)
res[ , -1] <- round(res[ , -1], 3)

cat("\n==== TEST-SET PERFORMANCE ====\n"); print(res, row.names = FALSE)
write.csv(res, file.path(PATHS$tables, "T4_test_performance.csv"), row.names = FALSE)

# RF importance (randomForest %IncMSE / IncNodePurity)
vi_rf <- as.data.frame(importance(rf))
vi_rf$var <- rownames(vi_rf)
vi_rf <- vi_rf[order(-vi_rf[["%IncMSE"]]), ]
write.csv(vi_rf, file.path(PATHS$tables, "T5_rf_importance.csv"), row.names = FALSE)
cat("\nTop 8 RF importance (%IncMSE):\n")
print(head(vi_rf[, c("var", "%IncMSE", "IncNodePurity")], 8), row.names = FALSE)

saveRDS(list(rf = rf, rf_ranger = rf_ranger, train = train, test = test,
             formula = fml, predictors = pv, best = best, best_ntree = best_ntree),
        file.path(PATHS$models, "rf_fit.rds"))
message("03_rf_train OK | ntree=", best_ntree, " | mtry=", best$mtry,
        " | min.node=", best$min.node.size)
