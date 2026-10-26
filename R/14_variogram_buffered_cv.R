# 14_variogram_buffered_cv.R - Variogram, buffered block CV, NNDM
# Covers the spatial autocorrelation behind the validation designs:
#   - estimate of the autocorrelation range (empirical variogram of RF residuals)
#   - block CV with buffer = range
#   - NNDM (Mila et al. 2022): leave-one-out that excludes the training points
#     closer than a distance threshold
#   - distribution of the sizes of the spatial blocks (K = 5, 15, 40)
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger); library(gstat); library(sp)})

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
d   <- dat$model_df
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv  <- predictor_vars(d)
f_rf <- as.formula(paste(TARGET,"~",paste(pv,collapse="+")))
f_sp <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),
                                          collapse="+")))

# ---- 1. Empirical variogram -------------------------------------------------
# Computed on the residuals of the NON-spatial RF (that is, after absorbing the
# covariates but before the coordinates). What is left is the residual spatial
# autocorrelation, which sizes the buffer of a buffered block CV.
set.seed(SEED)
rf0 <- ranger(f_rf, d, num.trees = 400, mtry = 6, min.node.size = 5,
              seed = SEED, num.threads = max(1L, parallel::detectCores()-1L))
res <- d[[TARGET]] - rf0$predictions   # OOB residuals (covariate-only RF)
sp_df <- data.frame(x = d$X_UTM33N, y = d$Y_UTM33N, r = res)
coordinates(sp_df) <- ~ x + y
# Variogram of the residuals, with fine lags and a checked fit:
#  - 250 m lags up to 10 km: the first falls at about 160 m and the range is
#    covered by several lags. Finer lags do not help: below about 100 m the
#    pairs are dominated by coincident points (350 points share their
#    position with another one) and by nearly coincident pairs, where the
#    semivariance of the OOB residuals exceeds the variance, and with 50-100 m
#    lags the range jumps to 3-8 km. Beyond 10 km the pairs link different
#    reliefs and the variogram is erratic;
#  - Sph, Exp and Gau from six starting values each: the fits that converge
#    without negative parameters are kept and the one with the smallest SSE of
#    the fit itself is chosen (attribute SSErr, np/h^2 weights as in gstat);
#  - practical range: Sph = a, Exp = 3a, Gau = sqrt(3) a;
#  - if no fit converges the script stops, rather than passing an
#    unidentified range downstream.
# The sensitivity to the lag width (150-300 m) goes to T29_variogram_lags.csv.
VG_WIDTH  <- 250
VG_CUTOFF <- 10000
VG_LAGS   <- c(150, 200, 250, 300)
s2_res    <- var(res)
practical_range <- function(vm) {
  st <- vm[vm$model != "Nug", ]
  switch(as.character(st$model), Exp = 3 * st$range,
         Gau = sqrt(3) * st$range, st$range)
}
fit_best <- function(v) {
  fits <- list()
  for (m in c("Sph", "Exp", "Gau")) for (r0 in c(300, 1000, 3000))
    for (nf in c(0.5, 0.9)) {
      w <- NULL
      f <- tryCatch(withCallingHandlers(
        fit.variogram(v, vgm(psill = (1 - nf) * s2_res, model = m, range = r0,
                             nugget = nf * s2_res), warn.if.neg = FALSE),
        warning = function(e) { w <<- conditionMessage(e)
                                invokeRestart("muffleWarning") }),
        error = function(e) NULL)
      ok <- !is.null(f) && is.null(w) && !isTRUE(attr(f, "singular")) &&
            all(f$psill >= 0) && all(f$range >= 0) && max(f$range) > 0
      if (ok) fits[[length(fits) + 1]] <- f
    }
  if (!length(fits)) return(NULL)
  fits[[which.min(sapply(fits, attr, "SSErr"))]]
}
v_emp <- variogram(r ~ 1, sp_df, cutoff = VG_CUTOFF, width = VG_WIDTH)
v_fit <- fit_best(v_emp)
if (is.null(v_fit))
  stop("14: no variogram fit converges with ", VG_WIDTH, " m lags")
best_name <- as.character(v_fit$model[v_fit$model != "Nug"])
range_m  <- practical_range(v_fit)
sill_m   <- sum(v_fit$psill)
nugget_m <- v_fit$psill[v_fit$model == "Nug"]
cat(sprintf(paste0("Variogramma (cov-only OOB resid, lag %d m fino a %d m, primo a %.0f m): ",
                   "model=%s  practical range=%.0f m  sill=%.0f  nugget=%.0f (%.0f%% of the sill)\n"),
            VG_WIDTH, VG_CUTOFF, v_emp$dist[1], best_name, range_m, sill_m,
            nugget_m, 100 * nugget_m / sill_m))
# Range with 3 decimals: 33_cv_geometry.R reads it back from here and
# multiplies it by 3, and with one decimal its effective range would differ by
# half a metre from the one used here for the buffers and for NNDM.
write.csv(data.frame(model = as.character(v_fit$model),
                     psill = round(v_fit$psill,3),
                     range = round(v_fit$range,3)),
          file.path(PATHS$tables, "T29_variogram_fit.csv"), row.names = FALSE)

# Sensitivity of the range to the lag width, same fit and same cutoff
lag_sens <- do.call(rbind, lapply(VG_LAGS, function(wd) {
  v <- variogram(r ~ 1, sp_df, cutoff = VG_CUTOFF, width = wd)
  f <- fit_best(v)
  data.frame(width_m = wd, n_lags = nrow(v), first_lag_m = round(v$dist[1]),
             model = if (is.null(f)) NA else as.character(f$model[f$model != "Nug"]),
             practical_range_m = if (is.null(f)) NA else round(practical_range(f)),
             nugget_share = if (is.null(f)) NA else
               round(f$psill[f$model == "Nug"] / sum(f$psill), 3))
}))
write.csv(lag_sens, file.path(PATHS$tables, "T29_variogram_lags.csv"),
          row.names = FALSE)
cat("\n==== Practical range by lag width (cutoff", VG_CUTOFF, "m) ====\n")
print(lag_sens, row.names = FALSE)

# Variogram plot: the residuals are those of RF, the model without
# coordinates (rf0 above).
png(file.path(PATHS$figs, "F50_variogram.png"), 1100, 700, res = 150)
plot(v_emp$dist, v_emp$gamma, pch = 19, col = "#2c7fb8",
     xlab = "Lag distance (m)", ylab = "Semivariance",
     main = sprintf("Empirical variogram of OOB residuals (RF), %d m lags\nfitted range = %.0f m",
                    VG_WIDTH, range_m))
abline(v = range_m, col = "red", lty = 2, lwd = 2)
g_dist <- seq(0, max(v_emp$dist), length.out = 200)
g_gamma <- variogramLine(v_fit, dist_vector = g_dist)$gamma
lines(g_dist, g_gamma, col = "#d73027", lwd = 2)
legend("bottomright", bty = "n",
       legend = c("empirical", sprintf("fitted (%s)", best_name),
                  sprintf("range = %.0f m", range_m)),
       col = c("#2c7fb8","#d73027","red"), pch = c(19, NA, NA),
       lty = c(NA, 1, 2), lwd = c(NA, 2, 2))
dev.off()

# ---- 2. Distribution of block sizes -----------------------------------------
block_sizes <- function(K) {
  set.seed(1)
  b <- kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]), K, nstart = 5,
              iter.max = 50)$cluster
  tab <- as.integer(table(b))
  data.frame(K = K, n_blocks = K, min = min(tab), q25 = quantile(tab, .25),
             median = median(tab), q75 = quantile(tab, .75), max = max(tab),
             cv = round(sd(tab)/mean(tab), 3))
}
bsz <- do.call(rbind, lapply(c(5, 15, 40), block_sizes))
write.csv(bsz, file.path(PATHS$tables, "T30_block_sizes.csv"),
          row.names = FALSE)
cat("\n==== Block-size distribution ====\n"); print(bsz, row.names = FALSE)

# ---- 3. Buffered block CV (K=15) --------------------------------------------
# For each test block: training = (all) \ (test block) \ (points within the
# buffer of the test block). Buffer = variogram range.
set.seed(1)
b15 <- kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]), 15, nstart = 5,
              iter.max = 50)$cluster
coords <- as.matrix(d[, c("X_UTM33N","Y_UTM33N")])

nn_dist <- function(idx_test, idx_train) {
  # nearest-neighbour distance from each test point to the closest training point
  if (!length(idx_train)) return(rep(Inf, length(idx_test)))
  ct <- coords[idx_test, , drop = FALSE]
  cr <- coords[idx_train, , drop = FALSE]
  apply(ct, 1, function(p) sqrt(min((cr[,1]-p[1])^2 + (cr[,2]-p[2])^2)))
}

buffered_cv <- function(buffer) {
  obs <- pred_rfsp <- c()
  excl <- integer(0)
  for (k in sort(unique(b15))) {
    teI <- which(b15 == k)
    # exclusions: training points within the buffer of some point of teI
    rest <- which(b15 != k)
    if (buffer > 0) {
      ct <- coords[teI, , drop = FALSE]
      keep_train <- rest[sapply(rest, function(i) {
        min((ct[,1] - coords[i,1])^2 + (ct[,2] - coords[i,2])^2) > buffer^2
      })]
    } else {
      keep_train <- rest
    }
    excl <- c(excl, length(rest) - length(keep_train))
    if (length(keep_train) < 100) next
    rf <- ranger(f_sp, d[keep_train,], num.trees = 400, mtry = 6,
                 min.node.size = 5, seed = SEED, num.threads = 4)
    p <- predict(rf, d[teI,])$predictions
    obs <- c(obs, d[teI, TARGET]); pred_rfsp <- c(pred_rfsp, p)
  }
  list(RMSE = sqrt(mean((obs - pred_rfsp)^2)),
       R2   = 1 - sum((obs - pred_rfsp)^2)/sum((obs - mean(obs))^2),
       n    = length(obs), mean_excl = mean(excl))
}

buf_levels <- c(0, round(range_m * 0.5), round(range_m), round(range_m * 2))
buf_res <- do.call(rbind, lapply(buf_levels, function(b) {
  r <- buffered_cv(b)
  data.frame(buffer_m = b, n_scored = r$n,
             mean_excluded = round(r$mean_excl, 1),
             RMSE = round(r$RMSE, 2), R2 = round(r$R2, 3))
}))
write.csv(buf_res, file.path(PATHS$tables, "T31_buffered_cv.csv"),
          row.names = FALSE)
cat("\n==== Buffered block CV (K=15) ====\n"); print(buf_res, row.names = FALSE)

# ---- 4. NNDM (Nearest Neighbour Distance Matching) --------------------------
# Simplified implementation: for each point, predict using a training set from
# which the points closer than `range_m` are removed. It is done in chunks of
# 1500 points for efficiency (random test subsample).
set.seed(SEED)
test_idx <- sample(seq_len(nrow(d)), min(1500, nrow(d)))
nndm_rmse <- numeric(0); nndm_obs <- nndm_pred <- c()
# group the test points in chunks of 100 to limit the number of fits
B <- ceiling(length(test_idx) / 100)
chunks <- split(test_idx, ceiling(seq_along(test_idx) / 100))
for (i in seq_along(chunks)) {
  ti <- chunks[[i]]
  ct <- coords[ti, , drop = FALSE]
  keep_train <- setdiff(seq_len(nrow(d)), ti)
  keep_train <- keep_train[sapply(keep_train, function(j) {
    min((ct[,1]-coords[j,1])^2 + (ct[,2]-coords[j,2])^2) > range_m^2
  })]
  if (length(keep_train) < 100) next
  rf <- ranger(f_sp, d[keep_train,], num.trees = 400, mtry = 6,
               min.node.size = 5, seed = SEED, num.threads = 4)
  p <- predict(rf, d[ti,])$predictions
  nndm_obs  <- c(nndm_obs,  d[ti, TARGET])
  nndm_pred <- c(nndm_pred, p)
}
nndm <- data.frame(method = "NNDM (range buffer)",
                    n_scored = length(nndm_obs),
                    RMSE = round(sqrt(mean((nndm_obs - nndm_pred)^2)), 2),
                    R2 = round(1 - sum((nndm_obs - nndm_pred)^2) /
                                sum((nndm_obs - mean(nndm_obs))^2), 3))
write.csv(nndm, file.path(PATHS$tables, "T32_nndm.csv"),
          row.names = FALSE)
cat("\n==== NNDM ====\n"); print(nndm, row.names = FALSE)

saveRDS(list(variogram = v_fit, range = range_m, lag_sensitivity = lag_sens,
             blocks = bsz, buffered = buf_res, nndm = nndm),
        file.path(PATHS$models, "buffered_cv.rds"))
message("14_variogram_buffered_cv OK | range=", round(range_m),
        " m | NNDM R2=", nndm$R2)
