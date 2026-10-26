# 10_operational_model.R - Operational model
# Thickness map as decision support: accurate where field data exist in the
# area, interpretable (e2tree on the covariate trend), with uncertainty and a
# flag for extrapolation. The heavy part, the scenario CV, is already in
# scenarios.rds (from 11_scenarios.R).
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger); library(gstat); library(sp)})

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
d   <- dat$model_df
grid <- dat$newdata
SC  <- readRDS(file.path(PATHS$models, "scenarios.rds"))
e2  <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))   # trend interpretation
# Covariates from predictor_vars() (00_setup.R), as in 01, 03 and 07: the 16
# variables of Table 2.
pv  <- predictor_vars(d)
theme_set(theme_minimal(base_size = 11))

getm <- function(sc, met, m) SC[[sc]][m, met]

# ---- 1. Curve of the local data required (central contribution) ---------------
# p = fraction of field data already present in the target area (probings)
curve <- data.frame(
  p = c(0, 0.10, 0.25, 0.50),
  scenario = c("sp15","oper15_p10","oper15_p25","oper15_p50"))
for (met in c("rf","rfsp","rfsp_rk")) {
  curve[[paste0("R2_", met)]]   <- sapply(curve$scenario, getm, met = met, m = "R2")
  curve[[paste0("RMSE_", met)]] <- sapply(curve$scenario, getm, met = met, m = "RMSE")
}
write.csv(curve, file.path(PATHS$tables, "T23_data_requirement_curve.csv"),
          row.names = FALSE)
cat("==== Accuracy vs % of local data ====\n"); print(round(curve[,-2],3))

# Figure F40 is drawn by R/13_figures_en.R; this script produces table T23.

# ---- 2. Validity domain (where the map is valid) -------------------------------
val <- data.frame(
  scenario = c("random5","sp40","sp15","sp5","leave_comune"),
  regime = c("(leakage: not an honest estimate)","fine interpolation",
             "medium interpolation","local extrapolation",
             "regional extrapolation"),
  R2_RFsp = sapply(c("random5","sp40","sp15","sp5","leave_comune"),
                   getm, met = "rfsp", m = "R2"),
  R2_GPT  = sapply(c("random5","sp40","sp15","sp5","leave_comune"),
                   getm, met = "gpt", m = "R2"))
write.csv(val, file.path(PATHS$tables, "T24_validity_domain.csv"), row.names = FALSE)
cat("\n==== Validity domain ====\n")
print(data.frame(scenario = val$scenario, R2_RFsp = round(val$R2_RFsp,3),
                  R2_GPT = round(val$R2_GPT,3)))

# Figure F41 is drawn by R/13_figures_en.R from the registry T47 of
# R/30_scenario_registry.R; this script produces table T24.

# ---- 3. Comparison with the empirical models on the same folds -----------------
bench <- data.frame(
  scenario = c("sp15","oper15_p25","random5"),
  RFsp_R2  = sapply(c("sp15","oper15_p25","random5"), getm, met="rfsp", m="R2"),
  RFsp_RMSE= sapply(c("sp15","oper15_p25","random5"), getm, met="rfsp", m="RMSE"),
  GPT_R2   = sapply(c("sp15","oper15_p25","random5"), getm, met="gpt", m="R2"),
  GPT_RMSE = sapply(c("sp15","oper15_p25","random5"), getm, met="gpt", m="RMSE"))
write.csv(bench, file.path(PATHS$tables, "T25_fair_benchmark.csv"), row.names = FALSE)
cat("\n==== Like-for-like benchmark, RFsp vs GPT (same folds) ====\n")
print(round(bench[,-1],2))

# ---- 4. Operational model and uncertainty (quantiles) --------------------------
# The delivered map covers the whole region and uses a forest without
# coordinates, with the physical covariates only (z0, morphometry, optical
# indices). The coordinates (RFsp) improve interpolation where field data
# exist in the area (curve of local data), but far from the data they tie the
# estimate to the nearest training leaves and detach it from the local
# covariates, inflating the values: for example unrealistic thicknesses in
# the Cilento, where z0 is low. With the covariates only, the values stay tied
# to z0 and to morphology and are sensible everywhere. RFsp is still fitted
# for the diagnostic scripts 25 and 26.
f_sp  <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),collapse="+")))
f_cov <- as.formula(paste(TARGET,"~",paste(pv,collapse="+")))
set.seed(SEED)
rfsp_all <- ranger(f_sp, d, num.trees = 500, mtry = 6, min.node.size = 5,
                   seed = SEED, quantreg = TRUE,
                   num.threads = max(1L, parallel::detectCores()-1L))
set.seed(SEED)
rfcov_all <- ranger(f_cov, d, num.trees = 500, mtry = 6, min.node.size = 5,
                    seed = SEED, quantreg = TRUE,
                    num.threads = max(1L, parallel::detectCores()-1L))

# Forest for the grid quantiles (q05, q50, q95, pi_width): more trees and
# larger leaves, so more samples per leaf and quantiles that are less coarse
# and more stable. Also with the covariates only, like the delivered mean.
set.seed(SEED)
rfcov_q <- ranger(f_cov, d, num.trees = 2000, mtry = 6, min.node.size = 20,
                  seed = SEED, quantreg = TRUE,
                  num.threads = max(1L, parallel::detectCores()-1L))

# Out-of-sample coverage of the 90% interval (15-block spatial CV), with the
# same configuration as the map (2000 trees, node size 20): we validate what we
# show. Comparing it with the 400/5 configuration shows that the coverage does
# not change appreciably; what changes is the coarseness of the quantiles.
set.seed(1); b15 <- kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]),15,nstart=5,
                           iter.max=50)$cluster
cover_cfg <- function(nt, ns) {
  cov <- c(); for (k in sort(unique(b15))) {
    m <- ranger(f_cov, d[b15!=k,], num.trees=nt, mtry=6, min.node.size=ns,
                quantreg=TRUE, seed=SEED, num.threads=4)
    q <- predict(m, d[b15==k,], type="quantiles",
                 quantiles=c(.05,.95))$predictions
    cov <- c(cov, mean(d[b15==k,TARGET] >= q[,1] & d[b15==k,TARGET] <= q[,2]))
  }
  mean(cov)
}
pi_cov_prev <- cover_cfg(400, 5)     # previous configuration
pi_cov      <- cover_cfg(2000, 20)   # mapped configuration (delivered)
cat(sprintf("\nOut-of-sample PI90%% coverage (spatial CV): %.3f (400/5) -> %.3f (2000/20)\n",
            pi_cov_prev, pi_cov))

# ---- 5. Grid: operational prediction, uncertainty and domain ------------------
stopifnot(all(pv %in% names(grid)))
G <- grid[stats::complete.cases(grid[, c(pv,"X_UTM33N","Y_UTM33N")]), ]
G$pred <- predict(rfcov_all, G)$predictions         # operational mean, covariates only
qg <- predict(rfcov_q, G, type="quantiles",         # quantiles from the larger forest
              quantiles=c(.05,.5,.95))$predictions
G$q05 <- qg[,1]; G$q50 <- qg[,2]; G$q95 <- qg[,3]
G$pi_width <- G$q95 - G$q05
# The map covers the whole region, without a mask: the covariates-only forest
# does not extrapolate through the coordinates and does not produce the
# artefacts that would call for one. in_domain stays as reliability
# information: it flags the few cells where a covariate, for example z0 on the
# volcanic edifices, falls outside the observed range. It is not used to mask
# the maps.
dom <- sapply(pv, function(v) range(d[[v]], na.rm=TRUE))
G$in_domain <- apply(sapply(pv, function(v)
  G[[v]] >= dom[1,v] & G[[v]] <= dom[2,v]), 1, all)
G$rule_zone <- factor(e2_assign(e2$tree, G)$node)   # trend interpretation

# Uncertainty field smoothed for display (pi_width only): it is an interval
# width, not a point estimate, and it does not enter the coverage computation,
# which uses the raw quantiles. Light Gaussian focal mean, the same as in
# 26_smooth_uncertainty.R.
suppressPackageStartupMessages(library(terra))
.r_pi  <- terra::rast(data.frame(x = G$X_UTM33N, y = G$Y_UTM33N,
                                 z = G$pi_width),   # full coverage
                      type = "xyz", crs = "EPSG:32633")
.gk    <- terra::focalMat(.r_pi, d = 1.5 * terra::res(.r_pi)[1], type = "Gauss")
.r_pis <- terra::mask(terra::focal(.r_pi, w = .gk, fun = "sum",
                                   na.rm = TRUE, na.policy = "omit"), .r_pi)
G$pi_width_smooth <- terra::extract(.r_pis, cbind(G$X_UTM33N, G$Y_UTM33N))[, 1]

cat(sprintf("Grid: %d points (full coverage), covariates within range %.1f%%\n",
            nrow(G), 100*mean(G$in_domain)))

# Maps F42-F45 are drawn by R/13_figures_en.R from the grid T27 produced here.

out <- G[,c("X_UTM33N","Y_UTM33N","pred","q05","q50","q95","pi_width",
            "pi_width_smooth","rule_zone","in_domain")]
write.csv(out, file.path(PATHS$tables,"T27_operational_grid.csv"),
          row.names=FALSE)
write.csv(data.frame(metric="PI90_coverage_spatialCV",value=round(pi_cov,3)),
          file.path(PATHS$tables,"T26_uncertainty.csv"), row.names=FALSE)
saveRDS(list(rfsp=rfsp_all, curve=curve, validity=val, bench=bench,
             pi_cov=pi_cov), file.path(PATHS$models,"operational.rds"))
message("10_operational_model OK | PI90 cov=", round(pi_cov,3),
        " | in-domain=", round(100*mean(G$in_domain),1), "%")
