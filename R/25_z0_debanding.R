# 25_z0_debanding.R - Diagnosis and reduction of the colour bands in the maps
# -----------------------------------------------------------------------------
# The causes, measured on the data:
#  1. z0 (ThickIso) draws the bands: in the training set it has only about 60
#     distinct integer values (discretised isopachs) and the forest splits the
#     space on thresholds of z0, so the band boundaries follow the isopachs;
#  2. the RF mean (pred) is instead almost continuous (about 12,000 distinct
#     values): for that panel the bands come mostly from the display (turbo
#     palette and 1000 m cells);
#  3. the quantiles really are stepped: q95 and pi_width have few distinct
#     values, as is normal for a quantile regression forest.
# Two experiments:
#  A. z0 smoothed after the fact with the model held fixed: almost no effect,
#     because the splits stay the same; the bands are therefore structural, not
#     a simple rounding effect;
#  B. model refitted with z0 made continuous, in the training set and on the
#     grid: measures how much the bands would fade if z0 were a continuous
#     field.
# Output: GeoTIFFs, a comparison figure, a table of the metrics.
source("R/00_setup.R")
suppressPackageStartupMessages({library(terra); library(ranger)})

CRS_UTM33N <- "EPSG:32633"
# The GeoTIFFs go to the output folder of the run (PTCAM_OUT), not to a fixed
# path.
outdir <- file.path(PATHS$out, "geotiff")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
Z0 <- "ThickIso"
PAL <- hcl.colors(100, "Viridis")

dat  <- readRDS(file.path(PATHS$models, "data_clean.rds"))
oper <- readRDS(file.path(PATHS$models, "operational.rds"))
rfsp <- oper$rfsp
pv   <- dat$predictors
d    <- dat$model_df
grid <- dat$newdata
G    <- grid[stats::complete.cases(grid[, c(pv, "X_UTM33N", "Y_UTM33N")]), ]
stopifnot(Z0 %in% names(G), Z0 %in% names(d))

cat(sprintf("z0 in the TRAINING set: %d distinct integer values over %d points\n",
            length(unique(d[[Z0]])), nrow(d)))

# ---- Continuous z0: a shared smooth isopach field (training + grid) ---------
# A z0 raster is built from the grid and smoothed with a Gaussian kernel; its
# continuous value is then resampled both at the grid points and at the
# training points.
rz0  <- terra::rast(data.frame(x = G$X_UTM33N, y = G$Y_UTM33N, z = G[[Z0]]),
                    type = "xyz", crs = CRS_UTM33N)
gk   <- terra::focalMat(rz0, d = 2 * terra::res(rz0)[1], type = "Gauss")
rz0s <- terra::focal(rz0, w = gk, fun = "sum", na.rm = TRUE, na.policy = "omit")
fillna <- function(v, orig) { v[is.na(v)] <- orig[is.na(v)]; v }
z0_grid_cont  <- fillna(terra::extract(rz0s, cbind(G$X_UTM33N, G$Y_UTM33N))[,1], G[[Z0]])
z0_train_cont <- fillna(terra::extract(rz0s, cbind(d$X_UTM33N, d$Y_UTM33N))[,1], d[[Z0]])
cat(sprintf("continuous z0: %d distinct values (grid)\n",
            length(unique(round(z0_grid_cont, 3)))))

# ---- Experiment A: fixed model, z0 smoothed at prediction time only ---------
G$predA_orig <- predict(rfsp, G)$predictions
GA <- G; GA[[Z0]] <- z0_grid_cont
G$predA_cont <- predict(rfsp, GA)$predictions

# ---- Experiment B: model refitted with continuous z0 ------------------------
f_sp <- as.formula(paste(TARGET, "~",
        paste(c(pv, "X_UTM33N", "Y_UTM33N"), collapse = "+")))
dC <- d; dC[[Z0]] <- z0_train_cont
set.seed(SEED)
rf_cont <- ranger(f_sp, dC, num.trees = 500, mtry = 6, min.node.size = 5,
                  seed = SEED, quantreg = TRUE,
                  num.threads = max(1L, parallel::detectCores() - 1L))
GB <- G; GB[[Z0]] <- z0_grid_cont
G$predB_cont <- predict(rf_cont, GB)$predictions
qB <- predict(rf_cont, GB, type = "quantiles", quantiles = c(.05,.95))$predictions
G$pi_widthB  <- qB[,2] - qB[,1]

# domain mask
dom <- sapply(pv, function(v) range(d[[v]], na.rm = TRUE))
G$in_domain <- apply(sapply(pv, function(v)
  G[[v]] >= dom[1,v] & G[[v]] <= dom[2,v]), 1, all)

# ---- Banding metrics --------------------------------------------------------
# n distinct values; plateau = % of 4-adjacent pairs with |d|<1 cm (flat
# interior); edge_jump = mean jump at the edges; roughness = mean 3x3 focal sd.
band_metrics <- function(col) {
  r <- terra::rast(data.frame(x = G$X_UTM33N, y = G$Y_UTM33N,
                              z = ifelse(G$in_domain, G[[col]], NA)),
                   type = "xyz", crs = CRS_UTM33N)
  M  <- terra::as.matrix(r, wide = TRUE)
  dd <- c(abs(M[,-1]-M[,-ncol(M)]), abs(M[-1,]-M[-nrow(M),])); dd <- dd[is.finite(dd)]
  fsd <- terra::focal(r, w = 3, fun = "sd", na.rm = TRUE)
  c(n_unique  = length(unique(round(G[[col]][G$in_domain], 4))),
    plateau   = mean(dd < 1),
    edge_jump = mean(dd[dd >= 1]),
    roughness = mean(terra::values(fsd), na.rm = TRUE))
}
res <- as.data.frame(rbind(
  predA_z0orig   = band_metrics("predA_orig"),   # baseline (= map of the paper)
  predA_z0smooth = band_metrics("predA_cont"),   # A: smoothing with the model held fixed
  predB_retrain  = band_metrics("predB_cont")))  # B: refitted with continuous z0

cat("\n==== BANDING: mean thickness ====\n"); print(round(res, 3))
cat(sprintf(
"\nExperiment A (fixed model): roughness %.2f -> %.2f cm  (%+.0f%%)  [expected ~zero: frozen splits]\n",
  res["predA_z0orig","roughness"], res["predA_z0smooth","roughness"],
  100*(res["predA_z0smooth","roughness"]/res["predA_z0orig","roughness"]-1)))
cat(sprintf(
"Experiment B (refitted): roughness %.2f -> %.2f cm  (%+.0f%%) ; plateau %.1f%% -> %.1f%%\n",
  res["predA_z0orig","roughness"], res["predB_retrain","roughness"],
  100*(res["predB_retrain","roughness"]/res["predA_z0orig","roughness"]-1),
  100*res["predA_z0orig","plateau"], 100*res["predB_retrain","plateau"]))

# ---- Output: GeoTIFFs + three-panel visual comparison -----------------------
mk_r <- function(col) terra::rast(
  data.frame(x = G$X_UTM33N, y = G$Y_UTM33N, z = ifelse(G$in_domain, G[[col]], NA)),
  type = "xyz", crs = CRS_UTM33N)
terra::writeRaster(mk_r("predA_orig"), file.path(outdir,"thickness_z0orig.tif"),
                   overwrite = TRUE, gdal = "COMPRESS=DEFLATE")
terra::writeRaster(mk_r("predB_cont"), file.path(outdir,"thickness_z0retrain.tif"),
                   overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

png(file.path(PATHS$figs, "F47_debanding_compare.png"),
    width = 1700, height = 620, res = 130)
op <- par(mfrow = c(1,3), mar = c(2,2,3,4))
rng <- range(c(G$predA_orig, G$predB_cont)[G$in_domain], na.rm = TRUE)
plot(mk_r("predA_orig"), main="(paper) z0 with 60 levels", range=rng, col=PAL)
plot(mk_r("predA_cont"), main="A) smoothing, fixed model", range=rng, col=PAL)
plot(mk_r("predB_cont"), main="B) refitted, continuous z0", range=rng, col=PAL)
par(op); dev.off()

write.csv(cbind(experiment = rownames(res), round(res,4)),
          file.path(PATHS$tables, "T28_debanding.csv"), row.names = FALSE)
cat(sprintf("\nOK -> %s/thickness_z0{orig,retrain}.tif , figs/F47_debanding_compare.png\n", outdir))
message("25_z0_debanding OK")
