# 03c_tradeoff_nloi.R - Frontier between forest accuracy and fidelity (nLoI)
# Fidelity is the structural measure of e2tree (nLoI), not a correlation
# between predictions. A very complex forest is accurate but hard to interpret
# (high nLoI); a regularised one is interpretable (low nLoI) and loses little
# accuracy.
source("R/00_setup.R")
suppressPackageStartupMessages(library(future)); plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
pv  <- fit$predictors; fml <- fit$formula
train <- fit$train[, c(TARGET, pv)]; test <- fit$test[, c(TARGET, pv)]
LEVEL <- 6; SUBN <- 2000

# Induction on the full training set: e2tree() needs the forest's training set
# as 'data' (see the note in R/03b_e2tree_setup.R). In every row of the
# frontier, forest, tree and nLoI are computed on the same training set.
S <- train

grid <- data.frame(
  config = c("RF predictive (nt1200,ns5)", "nt600 ns25", "nt400 ns50",
             "nt300 ns80 [selected]"),
  ntree  = c(1200, 600, 400, 300),
  ns     = c(5, 25, 50, 80))

out <- data.frame()
for (i in seq_len(nrow(grid))) {
  set.seed(SEED)
  # importance = TRUE as in 03b: permutation importance consumes random
  # numbers, so without this option, at equal seed, the forest would differ
  # from rf_e2 of 03b and nLoI would not match that of 04_e2tree.R for the
  # same configuration.
  rf <- randomForest(fml, train, ntree = grid$ntree[i],
                     nodesize = grid$ns[i], mtry = 10, importance = TRUE)
  acc_rmse <- rmse(test[[TARGET]], as.numeric(predict(rf, test)))
  acc_r2   <- r2(test[[TARGET]],  as.numeric(predict(rf, test)))
  nloi <- tryCatch({
    Dm <- createDisMatrix(rf, data = S, label = TARGET,
                          parallel = list(active = FALSE, no_cores = 1))
    t2 <- e2tree(fml, S, Dm, rf,
                 list(impTotal = 0.1, maxDec = 1e-6, n = 2, level = LEVEL))
    eValidation(S, t2, Dm, test = "measures", graph = FALSE,
                n_perm = 99, seed = SEED)$loi$nloi
  }, error = function(e) { message(grid$config[i], ": ",
                                   conditionMessage(e)); NA })
  out <- rbind(out, data.frame(config = grid$config[i],
    RF_RMSE = round(acc_rmse, 2), RF_R2 = round(acc_r2, 3),
    nLoI = round(nloi, 4)))
  cat(sprintf("%-26s | RF RMSE=%.1f R2=%.3f | nLoI=%.4f\n",
              grid$config[i], acc_rmse, acc_r2, nloi))
}
cat("\n(LOW nLoI = high structural fidelity between e2tree and the ensemble)\n")
print(out, row.names = FALSE)
write.csv(out, file.path(PATHS$tables, "T15_accuracy_fidelity_tradeoff.csv"),
          row.names = FALSE)
message("03c_tradeoff_nloi OK")
