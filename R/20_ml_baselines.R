# 20_ml_baselines.R - Modern baselines: XGBoost + spatial GAM
# Adds to RF/RFsp/RK and the three empirical models: XGBoost (with
# coordinates) and a GAM with a thin-plate spline on (x,y) + linear terms in
# the covariates. Same folds, for a like-for-like comparison.
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger); library(xgboost); library(mgcv)})

d <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv <- predictor_vars(d)
f_sp <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),
                                          collapse="+")))

# Covariate matrix for XGBoost (numeric)
mk_X <- function(df) as.matrix(df[, c(pv, "X_UTM33N","Y_UTM33N")])

# GAM formula: thin-plate spline on (x,y) + linear terms in the covariates
f_gam <- as.formula(paste(TARGET, "~ s(X_UTM33N, Y_UTM33N, k = 100) +",
                          paste(pv, collapse = " + ")))

blocks <- function(K){
  set.seed(1)
  kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]),K,nstart=5,iter.max=50)$cluster
}
set.seed(1); rnd <- sample(rep(1:5, length.out = nrow(d)))
muni <- as.integer(factor(d$localitiy))

eval_models <- function(folds, oper_p = 0) {
  obs <- pred_rfsp <- pred_xgb <- pred_gam <- c()
  for (k in sort(unique(folds))) {
    if (oper_p > 0) {
      ink <- which(folds == k); set.seed(k)
      add <- sample(ink, floor(oper_p * length(ink)))
      teI <- setdiff(ink, add)
      trI <- c(which(folds != k), add)
    } else {
      teI <- which(folds == k); trI <- which(folds != k)
    }
    if (length(teI) < 5) next
    tr <- d[trI, ]; te <- d[teI, ]
    # RFsp
    rf <- ranger(f_sp, tr, num.trees = 400, mtry = 6, min.node.size = 5,
                 seed = SEED, num.threads = 4)
    p_rfsp <- predict(rf, te)$predictions
    # XGBoost
    Xtr <- mk_X(tr); Xte <- mk_X(te)
    set.seed(SEED)
    xgb <- xgboost::xgb.train(
      data = xgboost::xgb.DMatrix(Xtr, label = tr[[TARGET]]),
      params = list(objective = "reg:squarederror",
                    eta = 0.05, max_depth = 6, subsample = 0.8,
                    colsample_bytree = 0.8, nthread = 4),
      nrounds = 600, verbose = 0)
    p_xgb <- predict(xgb, xgboost::xgb.DMatrix(Xte))
    # GAM with a spatial spline
    g <- tryCatch(mgcv::gam(f_gam, data = tr, method = "REML"),
                  error = function(e) NULL)
    p_gam <- if (!is.null(g)) predict(g, te) else rep(mean(tr[[TARGET]]), nrow(te))
    obs <- c(obs, te[[TARGET]])
    pred_rfsp <- c(pred_rfsp, p_rfsp)
    pred_xgb  <- c(pred_xgb, p_xgb)
    pred_gam  <- c(pred_gam, p_gam)
  }
  data.frame(
    model = c("RFsp","XGBoost (with coords)","GAM (TPS on coords)"),
    R2    = c(r2(obs, pred_rfsp), r2(obs, pred_xgb), r2(obs, pred_gam)),
    RMSE  = c(rmse(obs, pred_rfsp), rmse(obs, pred_xgb), rmse(obs, pred_gam)))
}

scenarios <- list(
  random5 = list(folds = rnd, p = 0),
  sp40    = list(folds = blocks(40), p = 0),
  sp15    = list(folds = blocks(15), p = 0),
  sp5     = list(folds = blocks(5),  p = 0),
  oper15_p25 = list(folds = blocks(15), p = 0.25),
  leave_municipality = list(folds = muni, p = 0))

all_b <- data.frame()
for (sc in names(scenarios)) {
  cat(sprintf("Scenario %s...\n", sc))
  r <- eval_models(scenarios[[sc]]$folds, scenarios[[sc]]$p)
  r$scenario <- sc
  r$R2 <- round(r$R2, 3); r$RMSE <- round(r$RMSE, 1)
  all_b <- rbind(all_b, r)
  print(r, row.names = FALSE)
}

write.csv(all_b, file.path(PATHS$tables, "T44_ml_baselines.csv"),
          row.names = FALSE)
saveRDS(all_b, file.path(PATHS$models, "ml_baselines.rds"))
message("20_ml_baselines OK | rows=", nrow(all_b))
