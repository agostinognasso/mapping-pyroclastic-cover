# 24_geotiff_export.R - Export the operational panels as georeferenced GeoTIFFs
# The PNG figures carry no spatial metadata. The georeferenced product is
# rebuilt from the source table, which is already in UTM 33N / WGS84
# (EPSG:32633), 1000 m cells. No warping, no loss.
source("R/00_setup.R")
suppressPackageStartupMessages(library(terra))

CRS_UTM33N <- "EPSG:32633"   # UTM zona 33N, WGS84 (Campania)
# The GeoTIFFs go to the output folder of the run (PTCAM_OUT), not to a fixed
# path.
outdir <- file.path(PATHS$out, "geotiff")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

G <- read.csv(file.path(PATHS$tables, "T27_operational_grid.csv"),
              check.names = FALSE)
G$in_domain <- as.logical(G$in_domain)

# Full coverage: RF with the covariates only does not extrapolate through the
# coordinates, so the layers cover the whole grid, without a mask. in_domain
# stays in the CSV as reliability information.
mask_out <- function(v) v

# Continuous fields: one GeoTIFF per panel
fields <- c(pred = "expected_thickness_cm",
            pi_width = "pi90_width_cm",
            q95 = "q95_thickness_cm",
            q05 = "q05_thickness_cm",
            q50 = "q50_thickness_cm")

for (col in names(fields)) {
  df <- data.frame(x = G$X_UTM33N, y = G$Y_UTM33N, z = mask_out(G[[col]]))
  r <- terra::rast(df, type = "xyz", crs = CRS_UTM33N)
  names(r) <- fields[[col]]
  f <- file.path(outdir, paste0(fields[[col]], ".tif"))
  terra::writeRaster(r, f, overwrite = TRUE,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  cat(sprintf("wrote %-26s  [%d x %d, res %.0f m]\n",
              basename(f), ncol(r), nrow(r), terra::res(r)[1]))
}

# Categorical field: e2tree rule zones (raster with an attribute table)
dfz <- data.frame(x = G$X_UTM33N, y = G$Y_UTM33N,
                  z = G$rule_zone)
rz <- terra::rast(dfz, type = "xyz", crs = CRS_UTM33N)
names(rz) <- "rule_zone"
levels(rz) <- data.frame(value = sort(unique(na.omit(dfz$z))),
                         rule = sort(unique(na.omit(dfz$z))))
fz <- file.path(outdir, "rule_zones.tif")
terra::writeRaster(rz, fz, overwrite = TRUE, datatype = "INT2S",
                   gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
cat(sprintf("wrote %-26s  [categorical, %d zones]\n",
            basename(fz), length(unique(na.omit(dfz$z)))))

# Multiband stack, handy for QGIS (a single file with all the continuous layers)
stk <- terra::rast(lapply(names(fields), function(col) {
  df <- data.frame(x = G$X_UTM33N, y = G$Y_UTM33N, z = mask_out(G[[col]]))
  r <- terra::rast(df, type = "xyz", crs = CRS_UTM33N); names(r) <- fields[[col]]; r
}))
terra::writeRaster(stk, file.path(outdir, "operational_stack.tif"),
                   overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))

cat(sprintf("\nOK -> %s\nCRS: %s | in-domain cells: %d/%d\n",
            outdir, CRS_UTM33N, sum(G$in_domain), nrow(G)))
message("24_geotiff_export OK")
