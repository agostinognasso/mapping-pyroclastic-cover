# 16_empirical_recal.R - Empirical models recalibrated per fold
# Comparing RFsp with GPT/SAPT/SEPT as published is formally correct but not
# very informative, because the empirical models have a scale bias (R^2 ~ -2.4
# for SAPT). Here a model z = a + b * f(x) is refitted on the empirical values
# for every training fold, and a + b * f is used on the test set. The
# comparison is then on functional form, not on a constant level.
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger)})

d <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv <- predictor_vars(d)
f_sp <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),
                                          collapse="+")))

# ---- CV designs used: random, sp40, sp15, sp5, oper15_p25, leave_mun -------
blocks <- function(K){
  set.seed(1)
  kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]), K, nstart=5, iter.max=50)$cluster
}
set.seed(1); rnd <- sample(rep(1:5, length.out = nrow(d)))
muni <- as.integer(factor(d$localitiy))

# ---- Function: for each fold, recalibrate the empirical model (a + b * f) ---
recal_empirical <- function(folds, oper_p = 0) {
  obs_all <- list()
  pred_emp_raw <- pred_emp_rec <- pred_rfsp <- list()
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
    pr_sp <- predict(rf, te)$predictions
    # Empirical models: as published + recalibrated
    emp_recal <- list()
    for (col in c("GPT_model","SEPT_model","SAPT_model")) {
      if (!col %in% names(tr)) next
      lm_e <- tryCatch(lm(tr[[TARGET]] ~ tr[[col]]),
                       error = function(e) NULL)
      if (is.null(lm_e) || any(is.na(coef(lm_e)))) {
        emp_recal[[col]] <- te[[col]]
      } else {
        emp_recal[[col]] <- coef(lm_e)[1] + coef(lm_e)[2] * te[[col]]
      }
    }
    obs_all[[length(obs_all)+1]] <- te[[TARGET]]
    pred_rfsp[[length(pred_rfsp)+1]] <- pr_sp
    for (col in c("GPT_model","SEPT_model","SAPT_model")) {
      if (!col %in% names(te)) next
      pred_emp_raw[[col]] <- c(pred_emp_raw[[col]], te[[col]])
      pred_emp_rec[[col]] <- c(pred_emp_rec[[col]], emp_recal[[col]])
    }
  }
  obs <- unlist(obs_all); pr_sp <- unlist(pred_rfsp)
  out <- data.frame(model = "RFsp", R2 = r2(obs, pr_sp),
                    RMSE = rmse(obs, pr_sp))
  for (col in c("GPT_model","SEPT_model","SAPT_model")) {
    if (!is.null(pred_emp_raw[[col]])) {
      out <- rbind(out, data.frame(
        model = paste0(col, "_raw"),
        R2 = r2(obs, pred_emp_raw[[col]]),
        RMSE = rmse(obs, pred_emp_raw[[col]])))
      out <- rbind(out, data.frame(
        model = paste0(col, "_recal"),
        R2 = r2(obs, pred_emp_rec[[col]]),
        RMSE = rmse(obs, pred_emp_rec[[col]])))
    }
  }
  out
}

scenarios <- list(
  random5 = list(folds = rnd, p = 0),
  sp40    = list(folds = blocks(40), p = 0),
  sp15    = list(folds = blocks(15), p = 0),
  sp5     = list(folds = blocks(5),  p = 0),
  oper15_p25 = list(folds = blocks(15), p = 0.25),
  leave_municipality = list(folds = muni, p = 0))

all_res <- data.frame()
for (sc in names(scenarios)) {
  r <- recal_empirical(scenarios[[sc]]$folds, scenarios[[sc]]$p)
  r$scenario <- sc
  all_res <- rbind(all_res, r)
  cat(sprintf("Scenario %s\n", sc))
  print(transform(r, R2 = round(R2,3), RMSE = round(RMSE,1)),
        row.names = FALSE)
}

all_res$R2 <- round(all_res$R2, 3)
all_res$RMSE <- round(all_res$RMSE, 1)
write.csv(all_res,
          file.path(PATHS$tables, "T35_empirical_recalibrated.csv"),
          row.names = FALSE)

cat("\n==== sp15: raw vs recal ====\n")
print(subset(all_res, scenario == "sp15"), row.names = FALSE)

saveRDS(all_res, file.path(PATHS$models, "empirical_recal.rds"))
message("16_empirical_recal OK | rows=", nrow(all_res))
