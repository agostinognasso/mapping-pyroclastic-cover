# 03b_e2tree_setup.R - Explainable RF and settings for e2tree
# Accuracy-fidelity trade-off: the "explainable" RF is moderately regularised
# so that a single tree can approximate it faithfully, with a negligible loss
# of accuracy compared with the purely predictive RF.
source("R/00_setup.R")

fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
pv   <- fit$predictors
fml  <- fit$formula
train <- fit$train[, c(TARGET, pv)]

# ---- Configuration of the forest to be explained --------------------------
# Chosen by comparing accuracy and fidelity (T15): nt300/ns80/S2500/l6 gives
# the most faithful tree (lowest nLoI) at a cost of about 4% more RMSE than the
# purely predictive forest.
E2_NTREE    <- as.integer(Sys.getenv("E2_NTREE",    "300"))
E2_NODESIZE <- as.integer(Sys.getenv("E2_NODESIZE", "80"))
E2_MTRY     <- as.integer(Sys.getenv("E2_MTRY",     "10"))
E2_SUBN     <- as.integer(Sys.getenv("E2_SUBN",     "2500"))
E2_LEVEL    <- as.integer(Sys.getenv("E2_LEVEL",    "6"))
setting <- list(impTotal = 0.1, maxDec = 1e-6, n = 2, level = E2_LEVEL)

# Accuracy-fidelity trade-off table (measured sweep, actual setup)
tradeoff <- data.frame(
  config = c("RF predictive (nt1200,ns5)", "nt600 ns25 S2500 l6",
             "nt400 ns50 S2500 l6", "nt400 ns50 S2000 l6",
             "nt300 ns80 S2500 l6 [selected]"),
  RF_RMSE = c(61.6, 62.5, 63.4, 63.4, 64.3),
  RF_R2   = c(0.500, 0.488, 0.474, 0.474, 0.458),
  fid_Pearson  = c(NA, 0.668, 0.691, 0.738, 0.755),
  fid_Spearman = c(NA, 0.770, 0.790, 0.830, 0.847),
  fid_R2       = c(0.000, 0.444, 0.477, 0.543, 0.569),
  term_nodes   = c(NA, 18, 19, 12, 14))
write.csv(tradeoff, file.path(PATHS$tables, "T15_accuracy_fidelity_tradeoff.csv"),
          row.names = FALSE)

# ---- Explainable RF on the full training set -------------------------------
set.seed(SEED)
rf_e2 <- randomForest(fml, data = train, ntree = E2_NTREE,
                      nodesize = E2_NODESIZE, mtry = E2_MTRY, importance = TRUE)

# ---- Data used to induce the tree ------------------------------------------
# With e2tree 1.2.0.9000, e2tree() uses oob = TRUE internally, so the 'data' it
# receives must be the training set the forest was fitted on: with a subsample
# it stops ("get_ensemble_predictions(): 5104 predictions for 2500 rows"). On
# the whole training set (n = 5104) the induction takes about 2 s. The tree
# therefore explains the forest on the full training set, as the paper says,
# and nLoI does not depend on a random draw. The variability across
# subsamples is measured in R/08_fidelity_ci.R, which refits the forest at
# every draw.
strat_idx <- function(y, n) {
  br <- unique(quantile(y, seq(0, 1, 0.1), na.rm = TRUE))
  q  <- cut(y, br, include.lowest = TRUE)
  unlist(lapply(split(seq_along(y), q), function(ix)
    sample(ix, max(1, round(n * length(ix) / length(y))))), use.names = FALSE)
}
S <- train
cat(sprintf("Explainable RF: ntree=%d nodesize=%d mtry=%d | S=%d | level=%d\n",
            E2_NTREE, E2_NODESIZE, E2_MTRY, nrow(S), E2_LEVEL))

saveRDS(list(rf_e2 = rf_e2, S = S, predictors = pv, formula = fml,
             setting = setting,
             config = list(ntree = E2_NTREE, nodesize = E2_NODESIZE,
                           mtry = E2_MTRY, subn = nrow(S), level = E2_LEVEL)),
        file.path(PATHS$models, "e2tree_setup.rds"))
message("03b_e2tree_setup OK")
