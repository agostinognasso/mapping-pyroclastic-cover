# 19_z0_ablation.R - Ablation of the covariate z0
# z0 (ThickIso) is itself an interpolation of isopachs and dominates the
# importance. The sensitivity of the results is quantified by removing it, on
# the same scenarios as 11_scenarios.R.
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger);library(gstat);library(sp)})

d <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv_full <- predictor_vars(d)
pv_no_z0 <- setdiff(pv_full, "ThickIso")
cat("Covariates, full =", length(pv_full),
    " | covariates without z0 =", length(pv_no_z0), "\n")

f_full <- as.formula(paste(TARGET, "~", paste(c(pv_full,"X_UTM33N","Y_UTM33N"),
                                              collapse="+")))
f_no_z0 <- as.formula(paste(TARGET, "~", paste(c(pv_no_z0,"X_UTM33N","Y_UTM33N"),
                                                collapse="+")))

blocks <- function(K){
  set.seed(1)
  kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]),K,nstart=5,iter.max=50)$cluster
}
muni <- as.integer(factor(d$localitiy))
set.seed(1); rnd <- sample(rep(1:5, length.out = nrow(d)))

cv_eval <- function(folds, fml, oper_p = 0) {
  obs_all <- pred_all <- c()
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
    rf <- ranger(fml, d[trI,], num.trees = 400, mtry = 6,
                 min.node.size = 5, seed = SEED, num.threads = 4)
    p  <- predict(rf, d[teI,])$predictions
    obs_all <- c(obs_all, d[teI, TARGET])
    pred_all <- c(pred_all, p)
  }
  c(R2 = r2(obs_all, pred_all), RMSE = rmse(obs_all, pred_all))
}

scenarios <- list(
  random5 = list(folds = rnd, p = 0),
  sp40    = list(folds = blocks(40), p = 0),
  sp15    = list(folds = blocks(15), p = 0),
  sp5     = list(folds = blocks(5),  p = 0),
  oper15_p25 = list(folds = blocks(15), p = 0.25),
  leave_municipality = list(folds = muni, p = 0))

out <- data.frame()
for (sc in names(scenarios)) {
  m_full  <- cv_eval(scenarios[[sc]]$folds, f_full,  scenarios[[sc]]$p)
  m_noz0  <- cv_eval(scenarios[[sc]]$folds, f_no_z0, scenarios[[sc]]$p)
  out <- rbind(out, data.frame(
    scenario = sc,
    R2_full = round(m_full["R2"], 3),
    R2_no_z0 = round(m_noz0["R2"], 3),
    delta_R2 = round(m_noz0["R2"] - m_full["R2"], 3),
    RMSE_full = round(m_full["RMSE"], 1),
    RMSE_no_z0 = round(m_noz0["RMSE"], 1)))
  cat(sprintf("%s | full R2=%.3f no_z0 R2=%.3f delta=%.3f\n",
              sc, m_full["R2"], m_noz0["R2"],
              m_noz0["R2"] - m_full["R2"]))
}

write.csv(out, file.path(PATHS$tables, "T43_z0_ablation.csv"),
          row.names = FALSE)
cat("\n==== z0 (ThickIso) ablation ====\n"); print(out, row.names = FALSE)

saveRDS(out, file.path(PATHS$models, "z0_ablation.rds"))
message("19_z0_ablation OK")
