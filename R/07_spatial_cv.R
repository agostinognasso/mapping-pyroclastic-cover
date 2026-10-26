# 07_spatial_cv.R - Spatial block CV (folds = disjoint groups of municipalities)
# Out-of-sample error of the forest without spatial leakage and, for each fold,
# the structural fidelity (nLoI) of e2tree.
source("R/00_setup.R")
suppressPackageStartupMessages(library(future)); plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

d  <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
su <- readRDS(file.path(PATHS$models, "e2tree_setup.rds"))$config
pv <- predictor_vars(d)
fml <- as.formula(paste(TARGET, "~", paste(pv, collapse = "+")))
nthr <- max(1L, parallel::detectCores() - 1L)

# --- Spatial folds: municipalities assigned to K blocks (balanced greedy) ------
K <- 5
loc_n <- sort(table(d$localitiy), decreasing = TRUE)
fold_of <- setNames(integer(length(loc_n)), names(loc_n))
load_k <- numeric(K)
for (L in names(loc_n)) { k <- which.min(load_k)
  fold_of[L] <- k; load_k[k] <- load_k[k] + loc_n[L] }
d$fold <- fold_of[d$localitiy]
cat("Points per spatial fold:\n"); print(table(d$fold))

strat_idx <- function(y, n) {
  br <- unique(quantile(y, seq(0, 1, 0.1), na.rm = TRUE))
  q  <- cut(y, br, include.lowest = TRUE)
  unlist(lapply(split(seq_along(y), q), function(ix)
    sample(ix, max(1, round(n * length(ix) / length(y))))), use.names = FALSE)
}

res <- data.frame(); fid <- data.frame()
for (k in 1:K) {
  tr <- d[d$fold != k, c(TARGET, pv)]
  te <- d[d$fold == k, c(TARGET, pv)]
  # Forest: spatial out-of-sample accuracy
  rf <- ranger(fml, tr, num.trees = su$ntree, min.node.size = su$nodesize,
               mtry = su$mtry, seed = SEED, num.threads = nthr)
  p  <- predict(rf, te)$predictions
  res <- rbind(res, data.frame(fold = k, n_test = nrow(te),
    RMSE = rmse(te[[TARGET]], p), MAE = mae(te[[TARGET]], p),
    R2 = r2(te[[TARGET]], p)))
  # e2tree: structural fidelity (nLoI) per fold, without tree predictions
  nloi <- tryCatch({
    set.seed(SEED)
    S <- tr[strat_idx(tr[[TARGET]], min(1500, nrow(tr))), ]
    rfo <- randomForest(fml, S, ntree = su$ntree, nodesize = su$nodesize,
                        mtry = su$mtry)
    Dm <- createDisMatrix(rfo, data = S, label = TARGET,
                          parallel = list(active = FALSE, no_cores = 1))
    t2 <- e2tree(fml, S, Dm, rfo,
                 list(impTotal = 0.1, maxDec = 1e-6, n = 2, level = su$level))
    eValidation(S, t2, Dm, test = "measures", graph = FALSE,
                n_perm = 99, seed = SEED)$loi$nloi
  }, error = function(e) { message("fold ", k, " nLoI: ",
                                   conditionMessage(e)); NA })
  fid <- rbind(fid, data.frame(fold = k, nLoI = nloi))
  cat(sprintf("fold %d | ENSEMBLE RMSE=%.1f MAE=%.1f R2=%.3f | e2tree nLoI=%.4f\n",
              k, res$RMSE[k], res$MAE[k], res$R2[k], fid$nLoI[k]))
}

agg <- data.frame(
  metric = c("RMSE", "MAE", "R2", "nLoI (e2tree fidelity)"),
  mean   = round(c(mean(res$RMSE), mean(res$MAE), mean(res$R2),
                   mean(fid$nLoI, na.rm = TRUE)), 3),
  sd     = round(c(sd(res$RMSE), sd(res$MAE), sd(res$R2),
                   sd(fid$nLoI, na.rm = TRUE)), 3))
cat("\n==== SPATIAL BLOCK CV (mean +/- sd over", K, "folds) ====\n")
print(agg, row.names = FALSE)
cat("For comparison: a random split gave RMSE~62-64 (optimistic because of spatial leakage)\n")
write.csv(res, file.path(PATHS$tables, "T17_spatialcv_folds.csv"), row.names = FALSE)
write.csv(agg, file.path(PATHS$tables, "T18_spatialcv_summary.csv"), row.names = FALSE)
message("07_spatial_cv OK | mean spatial RMSE=", round(mean(res$RMSE), 1))
