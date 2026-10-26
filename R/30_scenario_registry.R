# 30_scenario_registry.R - Single registry of the validation scenarios
#
# Every validation design has an acronym, a name and a description, and the
# LaTeX tables of the paper that report it are all generated from here. A
# scenario therefore has the same numbers and the same name everywhere:
# copying tables by hand from the CSV files of different runs leads to
# inconsistent values and names.
source("R/00_setup.R")
TB <- PATHS$tables
rd <- function(f) read.csv(file.path(TB, f), check.names = FALSE)
has <- function(f) file.exists(file.path(TB, f))

# ---- The registry -----------------------------------------------------------
# acr    : acronym used in all tables and figures
# key    : key under which the scenario is saved in the CSV files
# design : name of the validation design
# descr  : a one-line description
#
# Order of the designs, the same in all tables and figures: random, blocks
# from the coarsest to the finest, extrapolation designs, operational
# designs. The buffer of B15b is read from T31 (see below).
REG <- data.frame(stringsAsFactors = FALSE, rbind(
  c("R5",   "random5",      "random 5-fold",
    "random split of the measurement points"),
  c("B5",   "sp5",          "block CV, $K=5$",
    "coarse spatial folds; local extrapolation"),
  c("B15",  "sp15",         "block CV, $K=15$",
    "medium spatial folds; the scale at which the map is meant to be used"),
  c("B15b", "buffered545",  "buffered block, $K=15$, $b=BUFFER$~m",
    "as B15, with a buffer at the residual variogram range"),
  c("B40",  "sp40",         "block CV, $K=40$",
    "fine spatial folds, median $187$ measurement points each"),
  c("NNDM", "nndm",         "nearest-neighbour distance matching",
    "test--train distances matched to the prediction grid"),
  c("LOMO", "leave_comune", "leave-one-municipality-out",
    "one of the $15$ surveyed municipalities withheld at a time"),
  c("LOSO", "leave_sector", "leave-one-sector-out",
    "one contiguous regional sector withheld at a time"),
  c("OP10", "oper15_p10",   "B15 $+$ $10\\%$ in-area data",
    "$10\\%$ of each held-out spatial fold added to training"),
  c("OP25", "oper15_p25",   "B15 $+$ $25\\%$ in-area data",
    "$25\\%$ of each held-out spatial fold added to training"),
  c("OP50", "oper15_p50",   "B15 $+$ $50\\%$ in-area data",
    "$50\\%$ of each held-out spatial fold added to training")))
names(REG) <- c("acr", "key", "design", "descr")

# Classification into three regimes, with thresholds set in advance by this work.
regime_of <- function(r2) ifelse(is.na(r2), "--",
  ifelse(r2 >= 0.40, "green (operational use)",
  ifelse(r2 >= 0.20, "amber (in-area support)", "red (not usable)")))

# ---- Collect the values from a single run -----------------------------------
val <- rd("T24_validity_domain.csv")
REG$R2_RFsp <- val$R2_RFsp[match(REG$key, val$scenario)]
REG$R2_GPT  <- val$R2_GPT [match(REG$key, val$scenario)]

ben <- rd("T28_fair_benchmark_full.csv")
for (cl in c("RFsp_R2", "RFsp_RMSE", "GPT_R2", "SAPT_R2", "SEPT_R2"))
  REG[[cl]] <- ben[[cl]][match(REG$key, ben$scenario)]
# T28 and T24 must agree on the same scenario: this is the check that was
# missing. Any divergence beyond the tolerance is reported.
chk <- !is.na(REG$R2_RFsp) & !is.na(REG$RFsp_R2)
dv  <- abs(REG$R2_RFsp[chk] - REG$RFsp_R2[chk])
if (any(dv > 0.002))
  warning("REGISTRY: T24 and T28 diverge on ",
          paste(REG$acr[chk][dv > 0.002], collapse = ", "))
cat("Consistency check T24 vs T28: max |diff| =", round(max(dv, 0), 4), "\n")

# RFsp from T24 where available, otherwise from T28
REG$RFsp <- ifelse(is.na(REG$R2_RFsp), REG$RFsp_R2, REG$R2_RFsp)

# curve: RF (covariates only) for the scenarios that define it
cur <- rd("T23_data_requirement_curve.csv")
REG$RF <- cur$R2_rf[match(REG$key, cur$scenario)]
REG$RFsp_RK <- cur$R2_rfsp_rk[match(REG$key, cur$scenario)]
REG$RMSE_RF <- cur$RMSE_rf[match(REG$key, cur$scenario)]
# For the curve designs (B15, OP10, OP25, OP50) RFsp comes from T23 at full
# precision: the threshold must be applied to the unrounded value (OP10,
# 0.3997, is amber even though it prints as 0.400 at three decimals). The B
# designs come from T24, also at full precision.
iC <- match(REG$key, cur$scenario)
REG$RFsp <- ifelse(is.na(iC), REG$RFsp, cur$R2_rfsp[iC])

# buffered + NNDM
if (has("T31_buffered_cv.csv")) {
  bu <- rd("T31_buffered_cv.csv")
  # B15b is the third level of 14_variogram_buffered_cv.R: a buffer equal to
  # the practical range of the residual variogram, read from T31.
  i <- 3L
  REG$RFsp[REG$acr == "B15b"] <- bu$R2[i]
  # The range exceeds 1,000 m: thousands with a thin space, as in the rest of
  # the paper ($1\,206$).
  b_lab <- sub("^([0-9]+)([0-9]{3})$", "\\1\\\\,\\2",
               formatC(round(bu$buffer_m[i]), format = "d"))
  REG$design[REG$acr == "B15b"] <- sub("BUFFER", b_lab,
                                       REG$design[REG$acr == "B15b"],
                                       fixed = TRUE)
}
if (has("T32_nndm.csv")) REG$RFsp[REG$acr == "NNDM"] <- rd("T32_nndm.csv")$R2[1]
if (has("T34_leave_massif_cv.csv")) {
  lm <- rd("T34_leave_massif_cv.csv")
  p <- lm$R2[tolower(lm[[1]]) %in% c("pooled", "overall")]
  if (length(p) == 1) REG$RFsp[REG$acr == "LOSO"] <- p
}
# The RMSE of RFsp in cm next to the R^2, so the thresholds can also be read
# in centimetres (Section 4.1). From the same source as the R^2: T22 for the
# designs of 11_scenarios.R, T31, T32 and T34 for B15b, NNDM and LOSO.
t22 <- rd("T22_scenarios_RMSE.csv")
REG$RMSE_RFsp <- t22$rfsp[match(REG$key, t22$scenario)]
if (has("T31_buffered_cv.csv"))
  REG$RMSE_RFsp[REG$acr == "B15b"] <- rd("T31_buffered_cv.csv")$RMSE[3L]
if (has("T32_nndm.csv"))
  REG$RMSE_RFsp[REG$acr == "NNDM"] <- rd("T32_nndm.csv")$RMSE[1]
if (has("T34_leave_massif_cv.csv")) {
  lm <- rd("T34_leave_massif_cv.csv")
  p <- lm$RMSE[tolower(lm[[1]]) %in% c("pooled", "overall")]
  if (length(p) == 1) REG$RMSE_RFsp[REG$acr == "LOSO"] <- p
}
REG$regime <- regime_of(REG$RFsp)
REG$regime[REG$acr == "R5"] <- "not an operating regime (leakage reference)"

write.csv(REG, file.path(TB, "T47_scenario_registry.csv"), row.names = FALSE)

# ---- Emit the LaTeX tables -------------------------------------------------
f3 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3))
f1 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 1))

# "Operating regime" table: acronym + design + description + regime + R2
r <- REG[!is.na(REG$RFsp), ]
writeLines(c(
 "% tab_validity.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
 sprintf("%s & %s & %s & %s & %s & %s \\\\",
         r$acr, r$design, r$descr, r$regime, f3(r$RFsp), f1(r$RMSE_RFsp))),
 file.path(PATHS$paper, "tab_validity.tex"))

# ---- Table of the comparison on the same folds ------------------------------
# Columns: RFsp, OK, LM+RK, RF+RK, XGB, GAM, the three empirical models
# recalibrated per fold and the three raw empirical models. The recalibrated
# versions come from T35, the machine-learning models from T44.
b <- REG[!is.na(REG$RFsp_R2), ]
pick <- function(tab, scen, mod, col = "R2") {
  if (is.null(tab)) return(NA_real_)
  i <- which(tab$scenario == scen & tab$model == mod)
  if (length(i) != 1) NA_real_ else as.numeric(tab[[col]][i])
}
ALIAS <- c(leave_municipality = "leave_comune")
unalias <- function(x) ifelse(x %in% names(ALIAS), ALIAS[x], x)
rec <- if (has("T35_empirical_recalibrated.csv"))
         rd("T35_empirical_recalibrated.csv") else NULL
ml  <- if (has("T44_ml_baselines.csv")) rd("T44_ml_baselines.csv") else NULL
if (!is.null(ml) && !"model" %in% names(ml)) names(ml)[1] <- "model"
if (!is.null(ml))  ml$scenario  <- unalias(ml$scenario)
if (!is.null(rec)) rec$scenario <- unalias(rec$scenario)

ml_names <- c(XGB = "XGBoost (with coords)", GAM = "GAM (TPS on coords)")
for (m in names(ml_names)) b[[m]] <- vapply(b$key, function(k)
  pick(ml, k, ml_names[[m]]), numeric(1))
for (m in c("GPT", "SAPT", "SEPT")) {
  b[[paste0(m, "_rc")]]  <- vapply(b$key, function(k)
    pick(rec, k, paste0(m, "_model_recal")), numeric(1))
  b[[paste0(m, "_raw")]] <- vapply(b$key, function(k)
    pick(rec, k, paste0(m, "_model_raw")), numeric(1))
}
# fall back to the raw values of T28 where T35 does not cover the scenario
b$GPT_raw  <- ifelse(is.na(b$GPT_raw),  b$GPT_R2,  b$GPT_raw)
b$SAPT_raw <- ifelse(is.na(b$SAPT_raw), b$SAPT_R2, b$SAPT_raw)
b$SEPT_raw <- ifelse(is.na(b$SEPT_raw), b$SEPT_R2, b$SEPT_raw)

# Ordinary kriging: the standard geostatistical baseline, already computed by
# eval_fold() in 11_scenarios.R (column "ok" of T21). For a paper about spatial
# validation it is the first comparison a reader of this field expects.
sc <- rd("T21_scenarios_R2.csv")
b$OK <- sc$ok[match(b$key, sc$scenario)]
# RF+RK column: RF plus the kriging of its out-of-bag residuals (column "rk"
# of T21, same folds and same kriging as OK, from eval_fold() in
# 11_scenarios.R). It supports the statement of Section 5.3: whoever holds a
# local survey does better to interpolate it (OK) than to use it to correct
# the RF layer.
b$RF_RK <- sc$rk[match(b$key, sc$scenario)]
if (anyNA(b$RF_RK))
  stop("tab_bench: RF+RK missing in T21 for ",
       paste(b$acr[is.na(b$RF_RK)], collapse = ", "))
# LM+RK column: linear regression kriging (OLS regression on the 16
# covariates plus ordinary kriging of the residuals), same folds and same
# kriging as OK (37_regression_kriging.R, which checks that it reproduces OK
# of T21).
lk <- rd("T55_regression_kriging.csv")
b$LM_RK <- lk$R2_lm_rk[match(b$key, lk$scenario)]
if (anyNA(b$LM_RK))
  stop("tab_bench: LM+RK missing in T55 for ",
       paste(b$acr[is.na(b$LM_RK)], collapse = ", "), ": run 37")

writeLines(c(
 "% tab_bench.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
 # Raw SEPT included, as in Figure 6.
 sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
         b$acr, f3(b$RFsp_R2), f3(b$OK), f3(b$LM_RK), f3(b$RF_RK), f3(b$XGB),
         f3(b$GAM), f3(b$GPT_rc), f3(b$SAPT_rc), f3(b$SEPT_rc),
         f3(b$GPT_raw), f3(b$SAPT_raw), f3(b$SEPT_raw))),
 file.path(PATHS$paper, "tab_bench.tex"))

# Sampling-requirement curve table
c4 <- REG[!is.na(REG$RF), ]
c4$share <- c("0", "10", "25", "50")[match(c4$acr,
              c("B15", "OP10", "OP25", "OP50"))]
c4 <- c4[order(match(c4$acr, c("B15", "OP10", "OP25", "OP50"))), ]

# Absolute density and sampling geometry (R/31_curve_density.R). Without these
# two columns the curve is in relative units and cannot be read as guidance
# for a campaign: see the note in tab:curve.
if (has("T44_curve_density.csv")) {
  den <- rd("T44_curve_density.csv")
  i <- match(c(0, .10, .25, .50)[match(c4$acr, c("B15", "OP10", "OP25", "OP50"))],
             den$p)
  c4$dens <- den$dens_pooled[i]
  c4$nn   <- den$nn_median_m[i]
} else {
  c4$dens <- NA; c4$nn <- NA
}
f2 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2))
f0 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 0))
writeLines(c(
 "% tab_curve.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
 sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
         c4$share, c4$acr, f2(c4$dens), f0(c4$nn),
         f3(c4$RF), f3(c4$RFsp), f3(c4$RFsp_RK), f1(c4$RMSE_RF))),
 file.path(PATHS$paper, "tab_curve.tex"))

# Accuracy/fidelity frontier table (supplement). It used to be typed by hand in
# the supplement .tex, which is exactly how the divergences this registry
# exists to prevent arise.
if (has("T40_frontier_explicit.csv")) {
  fr <- rd("T40_frontier_explicit.csv")
  lbl <- c("accuracy-optimal (nt$1200$, ns$5$)", "nt$600$, ns$25$",
           "nt$400$, ns$50$", "nt$300$, ns$80$ (selected)")
  ci  <- ifelse(is.na(fr$nLoI_lo), "--",
                sprintf("$[%s,%s]$", f3(fr$nLoI_lo), f3(fr$nLoI_hi)))
  writeLines(c(
   "% tab_frontier.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
   sprintf("%s & %s & %s & %s & %s & %s \\\\",
           lbl[seq_len(nrow(fr))], f1(fr$RF_RMSE), f3(fr$RF_R2),
           f1(fr$RMSE_pct_loss), f3(fr$nLoI), ci)),
   file.path(PATHS$paper, "tab_frontier.tex"))
}

# Table of the summariser controls (CART surrogate + random partition)
if (has("T47_cart_surrogate.csv")) {
  cs <- rd("T47_cart_surrogate.csv")
  writeLines(c(
   "% tab_cart.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
   sprintf("%s & %d & %s & %s & %s \\\\",
           cs$summariser, cs$leaves, f3(cs$prox_corr), f3(cs$nLoI_bin),
           ifelse(is.na(cs$in_sample_R2), "--", f3(cs$in_sample_R2)))),
   file.path(PATHS$paper, "tab_cart.tex"))
}

# z0 ablation table
if (has("T43_z0_ablation.csv")) {
  z <- rd("T43_z0_ablation.csv")
  names(z)[1] <- "key"
  z$key <- unalias(z$key)
  z$acr <- REG$acr[match(z$key, REG$key)]
  z$acr[is.na(z$acr)] <- z$key[is.na(z$acr)]
  z <- z[order(match(z$acr, REG$acr)), ]
  cn <- names(z)
  full <- cn[grepl("full", cn)][1]; noz <- cn[grepl("no_?z0", cn)][1]
  dl   <- cn[grepl("^d|delta", cn)][1]
  writeLines(c(
   "% tab_z0.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
   sprintf("%s & %s & %s & %s \\\\", z$acr, f3(z[[full]]), f3(z[[noz]]),
           f3(z[[dl]]))),
   file.path(PATHS$paper, "tab_z0.tex"))
}

# Table of the buffered block CV, from the numbers of T31.
if (has("T31_buffered_cv.csv")) {
  bu  <- rd("T31_buffered_cv.csv")
  thn <- function(x) sub("^([0-9]+)([0-9]{3})$", "\\1\\\\,\\2",
                         formatC(round(x), format = "d"))
  writeLines(c(
   "% tab_buffered.tex -- generated by R/30_scenario_registry.R. Do not edit by hand.",
   sprintf("$%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\\\",
           thn(bu$buffer_m), thn(bu$n_scored), thn(bu$mean_excluded),
           f1(bu$RMSE), f3(bu$R2))),
   file.path(PATHS$paper, "tab_buffered.tex"))
}

cat("\n==== SCENARIO REGISTRY ====\n")
print(REG[, c("acr", "design", "regime", "RFsp", "RF")], row.names = FALSE)
message("30_scenario_registry OK -> tab_validity/bench/curve/z0/buffered.tex in ",
        PATHS$paper)
