# 33_cv_geometry.R - Geometry of the validation designs: how far from the
# training set a test point falls, and how many test points lie within the
# autocorrelation range.
#
# The paper states that a random split measures interpolation between nearly
# duplicate neighbours rather than the ability to predict in a new area. Here
# the statement becomes a measurement: for each CV design, the distance to the
# nearest training point and the fraction of test points within the effective
# range of the residual variogram. Within that radius the neighbour still
# carries information on the value to predict, and the CV measures local
# interpolation; beyond it, it measures extrapolation.
#
# Each design is replicated as in the scripts that use it:
#   random5      set.seed(1); sample(rep(1:5, ...))     rnd in 11_scenarios.R
#   sp40/15/5    set.seed(1); kmeans(scale(coord), K)   blocks() in 11_scenarios.R
#   leave_comune as.integer(factor(localitiy))          muni in 11_scenarios.R
#   leave_sector kmeans(M = 6) on municipal centroids   15_leave_massif.R
# For each fold k the training set is everything outside k and the test set is
# k. For each test point the Euclidean distance (UTM 33N, metres) to the
# nearest training point is computed, and the distances are pooled over all
# folds.
#
# The buffered designs (B15b, NNDM) are excluded: their geometry is fixed by
# construction, with no training point within the buffer.
#
# The effective range is read from T29_variogram_fit.csv with the rule of
# practical_range() in 14_variogram_buffered_cv.R (exponential 3 * range,
# Gaussian sqrt(3) * range).
#
# Output: T48_cv_geometry.csv
source("R/00_setup.R")

d  <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
CO <- c("X_UTM33N", "Y_UTM33N")
XY <- as.matrix(d[, CO])
n  <- nrow(d)

# ---- 1. Effective range of the variogram ------------------------------------
vg <- read.csv(file.path(PATHS$tables, "T29_variogram_fit.csv"))
st <- vg[vg$model != "Nug", ]
if (nrow(st) != 1L)
  stop("T29_variogram_fit.csv: expected a single structural model, found ",
       nrow(st))
# Same rule as practical_range() in 14: Exp = 3a, Gau = sqrt(3) a, Sph = a.
eff_range <- switch(as.character(st$model), Exp = 3 * st$range,
                    Gau = sqrt(3) * st$range, st$range)
message(sprintf("variogram: model %s, range %.1f m -> effective range %.1f m",
                st$model, st$range, eff_range))

# ---- 2. Designs -------------------------------------------------------------
blocks <- function(K) {
  set.seed(1)
  kmeans(scale(d[, CO]), K, nstart = 5, iter.max = 50)$cluster
}
set.seed(1)
rnd <- sample(rep(1:5, length.out = n))

# Regional sectors: k-means (M = 6) on the centroids of the surveyed
# municipalities, a replica of 15_leave_massif.R. There are 15 municipalities.
centro <- aggregate(cbind(X_UTM33N, Y_UTM33N) ~ localitiy, data = d, FUN = mean)
set.seed(SEED)
centro$massif <- kmeans(scale(centro[, CO]), centers = 6,
                        nstart = 25, iter.max = 100)$cluster
sector <- centro$massif[match(d$localitiy, centro$localitiy)]
message(sprintf("surveyed municipalities: %d | regional sectors: %d",
                nrow(centro), length(unique(sector))))

designs <- list(
  random5      = rnd,
  sp40         = blocks(40),
  sp15         = blocks(15),
  sp5          = blocks(5),
  leave_comune = as.integer(factor(d$localitiy)),
  leave_sector = sector)

# ---- 3. Distance to the nearest training point, pooled over all folds -------
nn_pooled <- function(folds) {
  unlist(lapply(sort(unique(folds)), function(k) {
    teI <- which(folds == k); trI <- which(folds != k)
    if (!length(teI) || !length(trI)) return(NULL)
    tr <- XY[trI, , drop = FALSE]; te <- XY[teI, , drop = FALSE]
    apply(te, 1, function(q)
      sqrt(min((tr[, 1] - q[1])^2 + (tr[, 2] - q[2])^2)))
  }), use.names = FALSE)
}

out <- do.call(rbind, lapply(names(designs), function(nm) {
  f <- designs[[nm]]
  dist <- nn_pooled(f)
  data.frame(
    design         = nm,
    n_folds        = length(unique(f)),
    n_test         = length(dist),
    nn_median_m    = round(median(dist), 1),
    nn_q1_m        = round(unname(quantile(dist, .25)), 1),
    nn_q3_m        = round(unname(quantile(dist, .75)), 1),
    pct_within_rng = round(100 * mean(dist <= eff_range), 1))
}))
out$eff_range_m <- round(eff_range, 1)

write.csv(out, file.path(PATHS$tables, "T48_cv_geometry.csv"), row.names = FALSE)

# ---- 4. LaTeX rows for the paper (acronyms and order of the registry) -------
# Order of the registry, from the coarsest design to the finest.
ACR <- c(random5 = "R5", sp5 = "B5", sp15 = "B15", sp40 = "B40",
         leave_comune = "LOMO", leave_sector = "LOSO")
# Thousands separated by the LaTeX thin space, as in the rest of the paper
# ($7\,293$). formatC(big.mark=) loses the backslash, so it is inserted by hand.
fmt <- function(x)
  sub("^([0-9]+)([0-9]{3})$", "\\1\\\\,\\2", formatC(round(x), format = "d"))
rows <- vapply(names(ACR), function(k) {
  r <- out[out$design == k, ]
  sprintf("%s & $%s$ & $%s$ & $%s$ & $%.1f$ \\\\",
          ACR[[k]], fmt(r$nn_median_m), fmt(r$nn_q1_m), fmt(r$nn_q3_m),
          r$pct_within_rng)
}, character(1))
writeLines(c("% tab_cvgeom.tex -- generated by R/33_cv_geometry.R. Do not edit by hand.",
             rows),
           file.path(PATHS$paper, "tab_cvgeom.tex"))

cat("\n#### Geometry of the validation designs ####\n")
cat(sprintf("effective variogram range: %.0f m | n = %d points\n\n",
            eff_range, n))
print(out, row.names = FALSE)
cat("\nReading: pct_within_rng is the percentage of test points that have at least one\n",
    "training point within the autocorrelation range. Close to 100 the design\n",
    "measures local interpolation, not transfer to a new area.\n",
    sep = "")

message("33_cv_geometry OK -> T48_cv_geometry.csv")
