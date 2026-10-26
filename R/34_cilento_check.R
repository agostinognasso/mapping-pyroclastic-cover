# 34_cilento_check.R - Numbers for the Cilento quoted in the text (Sections 4.7
# and 5.2)
#
# All on the same area, the ellipse CILENTO_ELLIPSE of 00_setup.R (the same as
# in Figure S10), and with the same statistic, the median over the grid cells.
# RF is the delivered model (pred and q95 of T27); RFsp is rfsp_all of
# 10_operational_model.R, fitted on all points with the coordinates.
#
# Output: T50_cilento_check.csv
source("R/00_setup.R")
suppressPackageStartupMessages(library(ranger))

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
opm <- readRDS(file.path(PATHS$models, "operational.rds"))
G   <- read.csv(file.path(PATHS$tables, "T27_operational_grid.csv"))
pv  <- predictor_vars(dat$model_df)

g <- dat$newdata[stats::complete.cases(
  dat$newdata[, c(pv, "X_UTM33N", "Y_UTM33N")]), ]
g <- merge(g, G[, c("X_UTM33N", "Y_UTM33N", "pred", "q95")],
           by = c("X_UTM33N", "Y_UTM33N"))
gc <- g[in_cilento(g$X_UTM33N, g$Y_UTM33N), ]
if (!nrow(gc)) stop("no grid cell in the Cilento ellipse")

sp_mean <- predict(opm$rfsp, gc)$predictions
sp_q95  <- predict(opm$rfsp, gc, type = "quantiles",
                   quantiles = 0.95)$predictions[, 1]
d <- dat$model_df
n_pts <- sum(in_cilento(d$X_UTM33N, d$Y_UTM33N))   # expected: no point

out <- data.frame(
  quantity = c("grid_cells", "measurement_points", "z0_median_cm",
               "RF_median_cm", "RF_q95_median_cm",
               "RFsp_median_cm", "RFsp_q95_median_cm"),
  value = c(nrow(gc), n_pts, round(median(gc$ThickIso), 1),
            round(median(gc$pred), 1), round(median(gc$q95), 1),
            round(median(sp_mean), 1), round(median(sp_q95), 1)))
write.csv(out, file.path(PATHS$tables, "T50_cilento_check.csv"),
          row.names = FALSE)
cat("\n==== Cilento (ellipse of Figure S10) ====\n")
print(out, row.names = FALSE)
message("34_cilento_check OK -> T50_cilento_check.csv")
