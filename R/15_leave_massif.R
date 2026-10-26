# 15_leave_massif.R - Leave-one-massif-out CV
# The "massifs" (regional sectors in the paper) are groups of municipalities
# obtained by clustering their UTM centroids. They replace
# leave-one-municipality-out (an administrative proxy) with a spatial scale
# consistent with the phenomenon (eruptive lobes / geological massifs).
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger); library(gstat); library(sp)})

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
d   <- dat$model_df
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv  <- predictor_vars(d)
f_sp <- as.formula(paste(TARGET,"~",paste(c(pv,"X_UTM33N","Y_UTM33N"),
                                          collapse="+")))

# ---- 1. Municipal centroids + clustering into M massifs ---------------------
centro <- aggregate(cbind(X_UTM33N, Y_UTM33N) ~ localitiy, data = d,
                    FUN = mean)
M <- 6
set.seed(SEED)
cl <- kmeans(scale(centro[, c("X_UTM33N","Y_UTM33N")]), centers = M,
             nstart = 25, iter.max = 100)$cluster
centro$massif <- cl
d$massif <- centro$massif[match(d$localitiy, centro$localitiy)]

cat("Points per massif:\n")
print(table(d$massif))

# Descriptive label (post hoc, from the barycentre)
mass_centro <- aggregate(cbind(X_UTM33N, Y_UTM33N) ~ massif, data = d,
                          FUN = mean)
mass_centro$label <- sprintf("M%d (x=%.0f, y=%.0f, n=%d)",
  mass_centro$massif, mass_centro$X_UTM33N, mass_centro$Y_UTM33N,
  as.integer(table(d$massif)[as.character(mass_centro$massif)]))
write.csv(mass_centro,
          file.path(PATHS$tables, "T33_massif_centroids.csv"),
          row.names = FALSE)
cat("\nMassifs (centroids):\n"); print(mass_centro, row.names = FALSE)

# Figure: points coloured by massif (ggplot, same style as the rest of the paper)
suppressPackageStartupMessages(library(ggplot2))
cols <- c("1"="#1b9e77","2"="#d95f02","3"="#7570b3",
          "4"="#e7298a","5"="#66a61e","6"="#e6ab02")
# Grey background with the study area, which also shows where data are
# missing (the central point of the paper); north arrow and scale bar; the
# volcanoes, Roccamonfina included; no title, which is in the caption.
grid_bg <- tryCatch({
  g <- read.csv(file.path(PATHS$tables, "T12_grid_predictions.csv"),
                check.names = FALSE)
  g[, c("X_UTM33N", "Y_UTM33N")]
}, error = function(e) NULL)

mg <- ggplot(d, aes(X_UTM33N, Y_UTM33N, colour = factor(massif)))
if (!is.null(grid_bg))
  mg <- mg + geom_raster(data = grid_bg, aes(X_UTM33N, Y_UTM33N),
                         inherit.aes = FALSE, fill = "grey88")
mg <- mg +
  geom_point(size = 0.5, alpha = 0.85) +
  geom_point(data = mass_centro, aes(X_UTM33N, Y_UTM33N), shape = 4,
             size = 4, stroke = 1.2, colour = "black", inherit.aes = FALSE) +
  geom_text(data = mass_centro,
            aes(X_UTM33N, Y_UTM33N, label = paste0("M", massif)),
            inherit.aes = FALSE, fontface = "bold", vjust = -0.9, size = 4) +
  scale_colour_manual(values = cols, guide = "none") +
  carto_layers(if (is.null(grid_bg)) d else grid_bg,
               vents = TRUE, scalebar_km = 20, north = TRUE,
               roccamonfina = TRUE, cilento = TRUE) +
  coord_equal() +
  labs(title = NULL, x = "X (UTM 33N, m)", y = "Y (UTM 33N, m)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
# Map: 600 dpi PNG (save_gg in 00_setup.R).
save_gg("F51_massif_groups.png", mg, 7.4, 6.2, vector = FALSE)

# ---- 2. Leave-one-massif-out CV ---------------------------------------------
n_thr <- max(1L, parallel::detectCores() - 1L)
results <- data.frame()
all_obs <- all_pred <- c()
for (m in sort(unique(d$massif))) {
  trI <- which(d$massif != m); teI <- which(d$massif == m)
  rf  <- ranger(f_sp, d[trI,], num.trees = 400, mtry = 6,
                min.node.size = 5, seed = SEED, num.threads = n_thr)
  p   <- predict(rf, d[teI,])$predictions
  obs <- d[teI, TARGET]
  results <- rbind(results, data.frame(
    massif = m, n_test = length(teI),
    RMSE = rmse(obs, p), MAE = mae(obs, p), R2 = r2(obs, p)))
  all_obs <- c(all_obs, obs); all_pred <- c(all_pred, p)
}

pooled <- data.frame(
  massif = "pooled", n_test = length(all_obs),
  RMSE = rmse(all_obs, all_pred), MAE = mae(all_obs, all_pred),
  R2 = r2(all_obs, all_pred))
results <- rbind(results, pooled)
results[, -1] <- lapply(results[, -1], function(x) {
  if (is.numeric(x)) round(x, 3) else x })
write.csv(results,
          file.path(PATHS$tables, "T34_leave_massif_cv.csv"),
          row.names = FALSE)
cat("\n==== Leave-one-massif-out CV ====\n")
print(results, row.names = FALSE)

saveRDS(list(results = results, massif = d$massif,
             centroids = mass_centro),
        file.path(PATHS$models, "leave_massif.rds"))
message("15_leave_massif OK | pooled R2=", pooled$R2)
