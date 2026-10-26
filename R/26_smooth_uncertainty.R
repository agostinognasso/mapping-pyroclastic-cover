# 26_smooth_uncertainty.R - A less blocky uncertainty panel
# -----------------------------------------------------------------------------
# The bands in the pi_width and q95 panels are real: a quantile regression
# forest returns few distinct values (about 912 for q95 and 1926 for
# pi_width, over 12,766 cells), because the quantiles are order statistics of
# the samples in the leaves. Two remedies, together:
#   (a) a richer quantile forest (more trees, larger min.node.size): more
#       samples per leaf, less coarse quantiles;
#   (b) a light focal mean on the uncertainty raster: the interval width is a
#       regular field and can be smoothed, point estimates cannot.
# The coverage does not change appreciably: the numerical coarseness is
# reduced, the level of the quantiles stays the same.
source("R/00_setup.R")
suppressPackageStartupMessages({library(terra); library(ranger)})

CRS_UTM33N <- "EPSG:32633"
# The GeoTIFFs go to the output folder of the run (PTCAM_OUT), not to a fixed
# path.
outdir <- file.path(PATHS$out, "geotiff")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat  <- readRDS(file.path(PATHS$models, "data_clean.rds"))
pv   <- dat$predictors
d    <- dat$model_df
grid <- dat$newdata
G    <- grid[stats::complete.cases(grid[, c(pv, "X_UTM33N", "Y_UTM33N")]), ]

f_sp <- as.formula(paste(TARGET, "~",
        paste(c(pv, "X_UTM33N", "Y_UTM33N"), collapse = "+")))

# ---- (a) Richer quantile forest ---------------------------------------------
set.seed(SEED)
qrf <- ranger(f_sp, d, num.trees = 2000, mtry = 6, min.node.size = 20,
              seed = SEED, quantreg = TRUE,
              num.threads = max(1L, parallel::detectCores() - 1L))
q  <- predict(qrf, G, type = "quantiles", quantiles = c(.05,.95))$predictions
G$pi_width_rich <- q[,2] - q[,1]

# baseline = panel of the paper (operational.rds)
oper <- readRDS(file.path(PATHS$models, "operational.rds"))
qb <- predict(oper$rfsp, G, type = "quantiles", quantiles = c(.05,.95))$predictions
G$pi_width_base <- qb[,2] - qb[,1]

# domain
dom <- sapply(pv, function(v) range(d[[v]], na.rm = TRUE))
G$in_domain <- apply(sapply(pv, function(v)
  G[[v]] >= dom[1,v] & G[[v]] <= dom[2,v]), 1, all)

mk_r <- function(col) terra::rast(
  data.frame(x = G$X_UTM33N, y = G$Y_UTM33N, z = ifelse(G$in_domain, G[[col]], NA)),
  type = "xyz", crs = CRS_UTM33N)

# ---- (b) Light focal smoothing of the uncertainty field ---------------------
r_base <- mk_r("pi_width_base")
r_rich <- mk_r("pi_width_rich")
gk <- terra::focalMat(r_rich, d = 1.5 * terra::res(r_rich)[1], type = "Gauss")
r_smooth <- terra::focal(r_rich, w = gk, fun = "sum", na.rm = TRUE, na.policy = "omit")
r_smooth <- terra::mask(r_smooth, r_rich)

# ---- Metrics: coarseness / roughness ----------------------------------------
band <- function(r) {
  M  <- terra::as.matrix(r, wide = TRUE)
  dd <- c(abs(M[,-1]-M[,-ncol(M)]), abs(M[-1,]-M[-nrow(M),])); dd <- dd[is.finite(dd)]
  v  <- terra::values(r); v <- v[is.finite(v)]
  fsd <- terra::focal(r, w = 3, fun = "sd", na.rm = TRUE)
  c(n_unique  = length(unique(round(v, 3))),
    plateau   = mean(dd < 1),
    roughness = mean(terra::values(fsd), na.rm = TRUE))
}
res <- as.data.frame(rbind(
  pi_base_paper   = band(r_base),    # 500 trees, node size 5 (paper)
  pi_rich_forest  = band(r_rich),    # 2000 trees, node size 20
  pi_rich_smooth  = band(r_smooth))) # + focal smoothing
cat("\n==== UNCERTAINTY (90% interval width): coarseness ====\n"); print(round(res,3))
cat(sprintf(
"\ndistinct values: %d -> %d -> %d   roughness: %.2f -> %.2f -> %.2f cm\n",
  res["pi_base_paper","n_unique"], res["pi_rich_forest","n_unique"],
  res["pi_rich_smooth","n_unique"],
  res["pi_base_paper","roughness"], res["pi_rich_forest","roughness"],
  res["pi_rich_smooth","roughness"]))

# ---- Output -----------------------------------------------------------------
terra::writeRaster(r_base,   file.path(outdir,"pi90_width_paper.tif"),
                   overwrite=TRUE, gdal="COMPRESS=DEFLATE")
terra::writeRaster(r_smooth, file.path(outdir,"pi90_width_smooth.tif"),
                   overwrite=TRUE, gdal="COMPRESS=DEFLATE")

suppressPackageStartupMessages(library(viridisLite))
png(file.path(PATHS$figs, "F49_uncertainty_smooth.png"),
    width = 1700, height = 620, res = 130)
op <- par(mfrow = c(1,3), mar = c(2,2,3,5))
rng <- range(c(terra::values(r_base), terra::values(r_smooth)), na.rm = TRUE)
pal <- rev(mako(100))
plot(r_base,   main="(paper) 500 trees, node size 5", range=rng, col=pal)
plot(r_rich,   main="(a) 2000 trees, node size 20",   range=rng, col=pal)
plot(r_smooth, main="(a)+(b) + focal smoothing", range=rng, col=pal)
par(op); dev.off()

write.csv(cbind(variant = rownames(res), round(res,4)),
          file.path(PATHS$tables, "T29_uncertainty_smooth.csv"), row.names = FALSE)
cat(sprintf("\nOK -> %s/pi90_width_{paper,smooth}.tif , figs/F49_uncertainty_smooth.png\n", outdir))
message("26_smooth_uncertainty OK")
