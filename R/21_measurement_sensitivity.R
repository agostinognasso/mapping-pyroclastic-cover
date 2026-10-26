# 21_measurement_sensitivity.R - Sensitivity to the measurement type
# PT-Cam-2 mixes probing tests ("fioretto"), boreholes ("sondaggio"), pits
# ("pozzetto"), trenches ("trincea"), outcrops, etc. (the labels in quotes are
# the values in the published data). This script checks that the results do
# not depend on the mixture, by restricting to the dominant type (probing
# tests, ~95%).
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger)})

d <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
cat("Measurement types:\n"); print(table(d$type))

# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv <- predictor_vars(d)
f_sp <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),
                                          collapse="+")))

blocks <- function(K, data){
  set.seed(1)
  kmeans(scale(data[,c("X_UTM33N","Y_UTM33N")]),K,nstart=5,iter.max=50)$cluster
}

cv_R2 <- function(data, folds, oper_p = 0) {
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
    rf <- ranger(f_sp, data[trI,], num.trees = 400, mtry = 6,
                 min.node.size = 5, seed = SEED, num.threads = 4)
    p <- predict(rf, data[teI,])$predictions
    obs_all <- c(obs_all, data[teI, TARGET])
    pred_all <- c(pred_all, p)
  }
  c(R2 = r2(obs_all, pred_all), RMSE = rmse(obs_all, pred_all))
}

# Subset: probing tests only ("fioretto", the dominant type)
d_fior <- d[d$type == "fioretto", ]
cat(sprintf("\nFull n=%d | fioretto-only n=%d\n", nrow(d), nrow(d_fior)))

scenarios <- list(
  sp15 = list(K = 15, p = 0),
  oper15_p25 = list(K = 15, p = 0.25))

out <- data.frame()
for (sc in names(scenarios)) {
  m_full  <- cv_R2(d,      blocks(scenarios[[sc]]$K, d),      scenarios[[sc]]$p)
  m_fior  <- cv_R2(d_fior, blocks(scenarios[[sc]]$K, d_fior), scenarios[[sc]]$p)
  out <- rbind(out, data.frame(
    scenario = sc,
    R2_all_types = round(m_full["R2"], 3),
    R2_fioretto_only = round(m_fior["R2"], 3),
    delta = round(m_fior["R2"] - m_full["R2"], 3),
    RMSE_all_types = round(m_full["RMSE"], 1),
    RMSE_fioretto_only = round(m_fior["RMSE"], 1)))
  cat(sprintf("%s | all R2=%.3f  fioretto-only R2=%.3f\n",
              sc, m_full["R2"], m_fior["R2"]))
}

write.csv(out, file.path(PATHS$tables, "T45_measurement_sensitivity.csv"),
          row.names = FALSE)
cat("\n==== Measurement-type sensitivity ====\n")
print(out, row.names = FALSE)
saveRDS(out, file.path(PATHS$models, "measurement_sensitivity.rds"))
message("21_measurement_sensitivity OK")
