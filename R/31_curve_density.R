# 31_curve_density.R - Absolute density (points/km^2) of the local-data curve.
#
# The curve of 11_scenarios.R is parametrised by p, the fraction of the points
# of the block moved to the training set (run_oper: add <- sample(ink,
# floor(p*length(ink)))). It is a relative quantity: to plan a campaign one
# needs to know how many points per km^2, and this script converts p into an
# absolute density.
#
# The blocks are the k-means on the scaled coordinates of the field points, as
# blocks() in 11_scenarios.R, with the same set.seed(1). The area of a block
# is obtained by assigning each cell of the prediction grid (1 km x 1 km) to
# the nearest block centroid, in the same scaled space as the k-means: this is
# the Voronoi partition induced by the k-means.
#
# Output: T44_curve_density.csv
source("R/00_setup.R")

d  <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
CO <- c("X_UTM33N", "Y_UTM33N")

# ---- 1. Blocks: exact replica of blocks(15) in 11_scenarios.R --------------
sc <- scale(d[, CO])                      # centre/scale to keep
set.seed(1)
km <- kmeans(sc, 15, nstart = 5, iter.max = 50)
b  <- km$cluster
stopifnot(length(b) == nrow(d))

# ---- 2. Block areas from the in-domain grid ---------------------------------
g <- read.csv(file.path(PATHS$tables, "T27_operational_grid.csv"))
g <- g[g$in_domain %in% c(TRUE, "TRUE"), ]

# grid resolution, measured (not assumed)
ux   <- sort(unique(g$X_UTM33N))
res  <- min(diff(ux))
cell_km2 <- (res / 1000)^2
message(sprintf("grid: resolution %.0f m -> %.2f km^2/cell | %d in-domain cells",
                res, cell_km2, nrow(g)))

# same affine transformation as applied to the field points
gs <- cbind((g$X_UTM33N - attr(sc, "scaled:center")[1]) / attr(sc, "scaled:scale")[1],
            (g$Y_UTM33N - attr(sc, "scaled:center")[2]) / attr(sc, "scaled:scale")[2])

# assignment to the nearest block centroid
d2 <- outer(gs[, 1], km$centers[, 1], "-")^2 +
      outer(gs[, 2], km$centers[, 2], "-")^2
g$block <- max.col(-d2, ties.method = "first")

area_km2 <- tabulate(g$block, nbins = 15) * cell_km2
n_pts    <- tabulate(b,       nbins = 15)

# ---- 2b. Area actually surveyed: convex hull of the points -------------------
# The Voronoi area is the territory of the sector; the convex hull of its field
# points is the part actually surveyed. The ratio between the two says how
# concentrated the points are, and is the reason why the density over the
# territory is not a recipe for a campaign.
hull_km2 <- sapply(1:15, function(k) {
  pk <- d[b == k, CO]
  if (nrow(pk) < 3) return(NA_real_)
  h <- chull(pk); x <- pk[h, 1]; y <- pk[h, 2]
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2 / 1e6   # shoelace, km^2
})

blk <- data.frame(block = 1:15, n_points = n_pts,
                  area_km2 = area_km2, hull_km2 = round(hull_km2, 1),
                  dens_territory = round(n_pts / area_km2, 3),
                  dens_surveyed  = round(n_pts / hull_km2, 3))
write.csv(blk, file.path(PATHS$tables, "T44_block_geometry.csv"), row.names = FALSE)

# ---- 3. In-area density for each p of the curve -----------------------------
# run_oper adds floor(p * n_k) points to block k: the count does not depend on
# the seed, only on p and on the number of points of the block.
P <- c(0, .10, .25, .50)
out <- do.call(rbind, lapply(P, function(p) {
  added <- floor(p * n_pts)
  dens  <- added / area_km2                  # density per block
  data.frame(
    p                = p,
    n_added_total    = sum(added),
    area_km2_total   = sum(area_km2),
    dens_pooled      = sum(added) / sum(area_km2),
    dens_median      = median(dens),
    dens_q1          = unname(quantile(dens, .25)),
    dens_q3          = unname(quantile(dens, .75)),
    # same density referred to the surveyed portion only
    dens_surveyed_pooled = sum(added) / sum(hull_km2, na.rm = TRUE))
}))

# ---- 4. At what distance does the model interpolate? -----------------------
# Needed to read the curve: the density over the territory says nothing about
# the geometry of the added points. The sampling of run_oper() is replicated
# and, for each test point, the distance to the nearest training point is
# measured.
nnd <- do.call(rbind, lapply(P, function(p) {
  dist_all <- unlist(lapply(sort(unique(b)), function(k) {
    ink <- which(b == k)
    set.seed(k)                                   # as in run_oper()
    add <- sample(ink, floor(p * length(ink)))
    teI <- setdiff(ink, add)
    if (length(teI) < 5) return(NULL)
    trI <- c(which(b != k), add)
    tr <- as.matrix(d[trI, CO]); te <- as.matrix(d[teI, CO])
    apply(te, 1, function(q) sqrt(min((tr[, 1] - q[1])^2 + (tr[, 2] - q[2])^2)))
  }))
  data.frame(p = p,
             nn_median_m = median(dist_all),
             nn_q1_m     = unname(quantile(dist_all, .25)),
             nn_q3_m     = unname(quantile(dist_all, .75)))
}))
out <- merge(out, nnd, by = "p")

out[] <- lapply(out, function(x) if (is.numeric(x)) round(x, 4) else x)
write.csv(out, file.path(PATHS$tables, "T44_curve_density.csv"), row.names = FALSE)

cat("\n#### Block geometry (K=15) ####\n")
print(blk, row.names = FALSE)
cat(sprintf("\ntotal in-domain area: %.0f km^2 | points: %d | mean density: %.3f pt/km^2\n",
            sum(area_km2), sum(n_pts), sum(n_pts) / sum(area_km2)))
cat("\n#### In-area density along the curve ####\n")
print(out, row.names = FALSE)

message("31_curve_density OK -> T44_curve_density.csv, T44_block_geometry.csv")
