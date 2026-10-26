# 38_applicability.R - Sampled zone and area of applicability of the map
# Two tools to say where the delivered map rests on the data.
#
# 1. Sampled zone: the cells whose distance from the nearest measurement point
#    is less than the practical range of the variogram of the RF residuals
#    (1,206 m, 14_variogram_buffered_cv.R). Outside this zone no training
#    point falls within the range: this is the condition of the B15b and NNDM
#    designs (R2 -0.03 and -0.28, red regime). The zone is drawn on Figure 9
#    (13_figures_en.R) from the polygon saved here.
# 2. Area of applicability (AOA; Meyer and Pebesma 2021, CAST package): the
#    "in-domain" grid of Table 7 checks the covariates one by one (each within
#    the range observed at the points), so it does not see combinations of
#    values never observed together. The AOA uses the dissimilarity index
#    (DI): the distance from the nearest training point in the space of the
#    standardised covariates weighted by the permutation importance of RF,
#    divided by the mean distance between training points. The threshold is
#    the upper whisker (Q3 + 1.5 IQR) of the DI of the training points
#    computed across the folds of the B15 spatial CV, as in Meyer and Pebesma
#    (2021).
# Outputs: T56 (areas), F59 (supplement), sampled_zone.rds and applicability.rds.
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger); library(sf); library(CAST)})
dcl  <- readRDS(file.path(PATHS$models, "data_clean.rds"))
d    <- dcl$model_df
pv   <- predictor_vars(d)
RANGE <- readRDS(file.path(PATHS$models, "buffered_cv.rds"))$range
G <- read.csv(file.path(PATHS$tables, "T27_operational_grid.csv"))
G$in_domain <- as.logical(G$in_domain)
nd <- dcl$newdata
G <- merge(G, nd[, c("X_UTM33N", "Y_UTM33N", pv)], by = c("X_UTM33N", "Y_UTM33N"))
stopifnot(nrow(G) == 16947, !anyNA(G[, pv]))

# ---- 1. Sampled zone ---------------------------------------------------------
pts <- st_as_sf(d, coords = c("X_UTM33N", "Y_UTM33N"), crs = 32633)
zone <- st_union(st_buffer(st_geometry(pts), RANGE, nQuadSegs = 12))
cells <- st_as_sf(G[, c("X_UTM33N", "Y_UTM33N")], coords = c("X_UTM33N", "Y_UTM33N"),
                  crs = 32633)
G$dist_nn <- as.numeric(st_distance(cells, pts[st_nearest_feature(cells, pts), ],
                                    by_element = TRUE))
G$sampled <- G$dist_nn <= RANGE
saveRDS(list(zone = zone, range = RANGE), file.path(PATHS$models, "sampled_zone.rds"))

# ---- 2. Area of applicability -------------------------------------------------
# Weights: permutation importance of RF in the configuration of the map (500
# trees, minimum node size 5, all points); negative importances count as 0.
rf_imp <- ranger(as.formula(paste(TARGET, "~", paste(pv, collapse = "+"))),
                 d[, c(TARGET, pv)], num.trees = 500, mtry = 6, min.node.size = 5,
                 importance = "permutation", seed = SEED,
                 num.threads = max(1L, parallel::detectCores() - 1L))
imp <- pmax(ranger::importance(rf_imp), 0)
W <- as.data.frame(t(imp[pv]))
set.seed(1)
b15 <- kmeans(scale(d[, c("X_UTM33N", "Y_UTM33N")]), 15, nstart = 5,
              iter.max = 50)$cluster
tdi <- trainDI(train = d[, pv], variables = pv, weight = W, CVtest = b15,
               verbose = FALSE)
ao  <- aoa(newdata = G[, pv], trainDI = tdi, verbose = FALSE)
G$DI  <- as.numeric(ao$DI)
G$AOA <- as.logical(ao$AOA)

# ---- Areas (cells of 1 km2) -----------------------------------------------------
A <- data.frame(
  what = c("grid", "in_domain", "aoa", "in_domain_and_aoa", "in_domain_not_aoa",
           "sampled_zone", "sampled_and_aoa", "aoa_threshold", "range_m"),
  value = c(nrow(G), sum(G$in_domain), sum(G$AOA), sum(G$in_domain & G$AOA),
            sum(G$in_domain & !G$AOA), sum(G$sampled), sum(G$sampled & G$AOA),
            round(tdi$threshold, 4), round(RANGE)))
print(A, row.names = FALSE)
write.csv(A, file.path(PATHS$tables, "T56_applicability.csv"), row.names = FALSE)
write.csv(data.frame(variable = pv, weight = round(as.numeric(imp[pv]), 2)),
          file.path(PATHS$tables, "T56_aoa_weights.csv"), row.names = FALSE)
saveRDS(list(grid = G[, c("X_UTM33N", "Y_UTM33N", "in_domain", "AOA", "DI",
                          "dist_nn", "sampled")],
             threshold = tdi$threshold, weights = imp),
        file.path(PATHS$models, "applicability.rds"))

# ---- F59: in-domain grid and AOA, in the supplement -----------------------------
lev <- c("inside", "outside")
DD <- rbind(
  data.frame(X_UTM33N = G$X_UTM33N, Y_UTM33N = G$Y_UTM33N,
             v = factor(ifelse(G$in_domain, "inside", "outside"), lev),
             panel = "(a) Covariate ranges (in-domain grid)"),
  data.frame(X_UTM33N = G$X_UTM33N, Y_UTM33N = G$Y_UTM33N,
             v = factor(ifelse(G$AOA, "inside", "outside"), lev),
             panel = "(b) Area of applicability"))
zone_sf <- st_sf(geometry = zone)
p59 <- ggplot(DD, aes(X_UTM33N, Y_UTM33N, fill = v)) +
  geom_raster() +
  geom_sf(data = zone_sf, inherit.aes = FALSE, fill = NA, colour = "black",
          linewidth = 0.25) +
  coord_sf(datum = NA, crs = 32633) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_manual(values = c(inside = "#9ecae1", outside = "#636363"),
                    name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "bottom",
        strip.text = element_text(hjust = 0, size = 10)) +
  carto_layers(DD, text_size = 2.8, roccamonfina = TRUE, cilento = TRUE)
save_gg("F59_applicability.png", p59, 11, 5.6, vector = FALSE)
message("38_applicability OK")
