# 17_pi_calibration.R - Calibration of the prediction intervals (Section 3.4)
# Marginal 90% coverage is not enough: coverage within classes (by quartile of
# the prediction, slope class and sector) and a proper score such as the
# Winkler interval score are also needed.
#
# What is delivered is what is validated: the intervals of Figure 9 and of the
# GeoTIFFs come from RF with the covariates only and the quantile
# configuration of the map (rfcov_q of 10_operational_model.R: 2000 trees, node
# size 20), so T36-T39, Table 9 and Figure 7 concern RF. RFsp stays as a
# comparison in the supplement (T49, tab_supp_pi.tex), on the same 15 folds.
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger)})

d <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv <- predictor_vars(d)
f_cov <- as.formula(paste(TARGET,"~",paste(pv, collapse="+")))
f_sp  <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),
                                           collapse="+")))
n_thr <- max(1L, parallel::detectCores() - 1L)

# Configurations: delivered RF (as rfcov_q) and RFsp as the comparison
CFG <- list(
  RF   = list(f = f_cov, nt = 2000, ns = 20, label = "RF (delivered)"),
  RFsp = list(f = f_sp,  nt = 400,  ns = 5,  label = "RFsp (benchmark)"))

# Massif of each point (from 15_leave_massif.R)
if (!file.exists(file.path(PATHS$models, "leave_massif.rds")))
  stop("Run 15_leave_massif.R first")
lm_fit <- readRDS(file.path(PATHS$models, "leave_massif.rds"))
d$massif <- lm_fit$massif

# ---- Spatial CV with K=15 + quantiles 0.05/0.5/0.95 -------------------------
set.seed(1)
b15 <- kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]), 15, nstart = 5,
              iter.max = 50)$cluster

alpha <- 0.10
winkler <- function(y, lo, hi, a) {
  w <- hi - lo
  w + ifelse(y < lo, 2/a * (lo - y), 0) +
      ifelse(y > hi, 2/a * (y - hi), 0)
}
pinball <- function(y, q, p) ifelse(y >= q, (y-q)*p, (q-y)*(1-p))

run_cfg <- function(cf) {
  o <- list(obs = c(), q05 = c(), q50 = c(), q95 = c(), pred = c(),
            slope = c(), massif = c())
  for (k in sort(unique(b15))) {
    trI <- which(b15 != k); teI <- which(b15 == k)
    m <- ranger(cf$f, d[trI,], num.trees = cf$nt, mtry = 6,
                min.node.size = cf$ns, quantreg = TRUE, seed = SEED,
                num.threads = n_thr)
    q <- predict(m, d[teI,], type = "quantiles",
                 quantiles = c(0.05, 0.5, 0.95))$predictions
    o$obs <- c(o$obs, d[teI, TARGET])
    o$q05 <- c(o$q05, q[,1]); o$q50 <- c(o$q50, q[,2]); o$q95 <- c(o$q95, q[,3])
    o$pred   <- c(o$pred, predict(m, d[teI,])$predictions)
    o$slope  <- c(o$slope, d[teI, "Slope"])
    o$massif <- c(o$massif, d[teI, "massif"])
  }
  o$inside <- o$obs >= o$q05 & o$obs <= o$q95
  o$ws     <- winkler(o$obs, o$q05, o$q95, alpha)
  o$width  <- o$q95 - o$q05
  o
}
R <- lapply(CFG, run_cfg)

# ---- Marginal coverage + Winkler score -----------------------------------
marg_of <- function(o)
  data.frame(metric = c("coverage_90", "winkler_90", "mean_width_cm",
                        "pinball_q05", "pinball_q50", "pinball_q95"),
             value = round(c(mean(o$inside), mean(o$ws), mean(o$width),
                             mean(pinball(o$obs, o$q05, 0.05)),
                             mean(pinball(o$obs, o$q50, 0.50)),
                             mean(pinball(o$obs, o$q95, 0.95))), 3))
by_strata <- function(o, g, extra = FALSE) {
  out <- data.frame(stratum = levels(g), n = as.integer(table(g)))
  if (extra) out$mean_yhat <- round(tapply(o$pred, g, mean), 1)
  out$coverage_90 <- round(tapply(o$inside, g, mean), 3)
  out$winkler_90  <- round(tapply(o$ws, g, mean), 1)
  out$width       <- round(tapply(o$width, g, mean), 1)
  out
}
qh_of <- function(o) cut(o$pred, quantile(o$pred, c(0, .25, .5, .75, 1)),
                         include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))
sl_of <- function(o) cut(o$slope, c(-Inf, 5, 15, 30, Inf),
                         labels = c("flat <=5","gentle 5-15","steep 15-30",
                                    "very steep >30"))
ms_of <- function(o) factor(paste0("M", o$massif),
                            levels = paste0("M", sort(unique(o$massif))))

S <- lapply(R, function(o) list(
  marg = marg_of(o),
  by_yhat = by_strata(o, qh_of(o), extra = TRUE),
  by_slope = by_strata(o, sl_of(o)),
  by_massif = by_strata(o, ms_of(o))))

# T36-T39: delivered model
rf <- S$RF
cat("==== Marginal calibration (RF, delivered) ====\n")
print(rf$marg, row.names = FALSE)
write.csv(rf$marg, file.path(PATHS$tables, "T36_pi_marginal.csv"), row.names = FALSE)
cat("\n==== Coverage by yhat quartile (RF) ====\n"); print(rf$by_yhat, row.names = FALSE)
write.csv(rf$by_yhat, file.path(PATHS$tables, "T37_pi_by_yhat.csv"), row.names = FALSE)
cat("\n==== Coverage by slope (RF) ====\n"); print(rf$by_slope, row.names = FALSE)
write.csv(rf$by_slope, file.path(PATHS$tables, "T38_pi_by_slope.csv"), row.names = FALSE)
cat("\n==== Coverage by regional sector (RF) ====\n"); print(rf$by_massif, row.names = FALSE)
write.csv(rf$by_massif, file.path(PATHS$tables, "T39_pi_by_massif.csv"), row.names = FALSE)

# T49: RF / RFsp comparison on the same folds (supplement)
cmp <- do.call(rbind, lapply(names(S), function(m) rbind(
  data.frame(model = m, block = "marginal", stratum = S[[m]]$marg$metric,
             n = NA, value = S[[m]]$marg$value),
  data.frame(model = m, block = "regional_sector_coverage",
             stratum = S[[m]]$by_massif$stratum, n = S[[m]]$by_massif$n,
             value = S[[m]]$by_massif$coverage_90))))
write.csv(cmp, file.path(PATHS$tables, "T49_pi_rf_vs_rfsp.csv"), row.names = FALSE)
cat("\n==== RF vs RFsp (same folds) ====\n"); print(cmp, row.names = FALSE)

# ---- Figure 7 ---------------------------------------------------------------
# Figure 7 (F52_conditional_coverage) is drawn by 39_coverage_figure.R from
# pi_calibration.rds, with the Wilson intervals and the number of points of
# each class.

# ---- LaTeX rows: Table 9 (main text) and the comparison in the supplement ---
f3 <- function(x) formatC(x, format = "f", digits = 3)
f1 <- function(x) formatC(x, format = "f", digits = 1)
mv <- function(m, k) S[[m]]$marg$value[S[[m]]$marg$metric == k]
writeLines(c(
  "% tab_unc.tex -- generated by R/17_pi_calibration.R. Do not edit by hand.",
  sprintf("Coverage of the $90\\%%$ interval & $%s$\\\\", f3(mv("RF", "coverage_90"))),
  sprintf("Winkler interval score (cm) & $%s$\\\\", f1(mv("RF", "winkler_90"))),
  sprintf("Mean interval width (cm) & $%s$\\\\", f1(mv("RF", "mean_width_cm"))),
  sprintf("Pinball loss, $q_{0.05}$ & $%s$\\\\", f1(mv("RF", "pinball_q05"))),
  sprintf("Pinball loss, $q_{0.50}$ & $%s$\\\\", f1(mv("RF", "pinball_q50"))),
  sprintf("Pinball loss, $q_{0.95}$ & $%s$\\\\", f1(mv("RF", "pinball_q95")))),
  file.path(PATHS$paper, "tab_unc.tex"))

thn <- function(x) sub("^([0-9]+)([0-9]{3})$", "\\1\\\\,\\2",
                       formatC(round(x), format = "d"))   # 4722 -> 4\,722
sec <- S$RF$by_massif$stratum
cov_s <- function(m) S[[m]]$by_massif$coverage_90[match(sec, S[[m]]$by_massif$stratum)]
writeLines(c(
  "% tab_supp_pi.tex -- generated by R/17_pi_calibration.R. Do not edit by hand.",
  sprintf("Coverage of the $90\\%%$ interval & -- & $%s$ & $%s$\\\\",
          f3(mv("RF", "coverage_90")), f3(mv("RFsp", "coverage_90"))),
  sprintf("Winkler interval score (cm) & -- & $%s$ & $%s$\\\\",
          f1(mv("RF", "winkler_90")), f1(mv("RFsp", "winkler_90"))),
  sprintf("Mean interval width (cm) & -- & $%s$ & $%s$\\\\",
          f1(mv("RF", "mean_width_cm")), f1(mv("RFsp", "mean_width_cm"))),
  "\\midrule",
  sprintf("Coverage, regional sector %s & $%s$ & $%s$ & $%s$\\\\",
          sec, thn(S$RF$by_massif$n), f3(cov_s("RF")), f3(cov_s("RFsp")))),
  file.path(PATHS$paper, "tab_supp_pi.tex"))

saveRDS(S, file.path(PATHS$models, "pi_calibration.rds"))
message("17_pi_calibration OK | RF cov=", mv("RF", "coverage_90"),
        " Winkler=", mv("RF", "winkler_90"),
        " | RFsp cov=", mv("RFsp", "coverage_90"))
