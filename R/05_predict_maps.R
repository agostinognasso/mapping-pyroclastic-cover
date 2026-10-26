# 05_predict_maps.R - Out-of-sample maps for operational use
# The prediction on the grid always comes from the forest. e2tree does not
# predict: it is overlaid only as an interpretation, showing which rule each
# cell belongs to.
source("R/00_setup.R")

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
e2  <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))
rf  <- e2$rf            # the interpreted forest, also used to predict
pv  <- fit$predictors
grid <- dat$newdata

stopifnot(all(pv %in% names(grid)))
G <- grid[stats::complete.cases(grid[, pv]), ]
cat("Valid grid points:", nrow(G), "/", nrow(grid), "\n")

# ---- Prediction: the forest only ---------------------------------------------
G$pred_rf <- as.numeric(predict(rf, G))

# ---- Extrapolation mask (outside the training domain) -------------------------
trn_ref <- fit$train
dom <- sapply(pv, function(v) range(trn_ref[[v]], na.rm = TRUE))
G$in_domain <- apply(sapply(pv, function(v)
  G[[v]] >= dom[1, v] & G[[v]] <= dom[2, v]), 1, all)
G$extrap_vars <- apply(sapply(pv, function(v)
  G[[v]] < dom[1, v] | G[[v]] > dom[2, v]), 1, sum)
cat(sprintf("In-domain points: %d / %d (%.1f%%)\n",
            sum(G$in_domain), nrow(G), 100 * mean(G$in_domain)))

# ---- Interpretation with e2tree: the rule of each cell -------------------------
# Only which rule applies, as a label: no numeric value estimated by the tree,
# which is not a predictor.
zn <- e2_assign(e2$tree, G)
G$rule_zone <- factor(zn$node)

cat("RF (ensemble) grid:",
    paste(round(summary(G$pred_rf), 1), collapse = " "), "\n")

# ---- Thickness classes of the forest (Jenks/Fisher, as in the paper) ----------
if (!requireNamespace("classInt", quietly = TRUE))
  try(install.packages("classInt", repos = "https://cloud.r-project.org",
                        quiet = TRUE), silent = TRUE)
jenks_breaks <- function(v, k = 5) {
  if (requireNamespace("classInt", quietly = TRUE))
    unique(classInt::classIntervals(v, n = k, style = "fisher")$brks)
  else unique(stats::quantile(v, seq(0, 1, length.out = k + 1), na.rm = TRUE))
}
br <- jenks_breaks(G$pred_rf, 5)
G$class_rf <- cut(G$pred_rf, breaks = br, include.lowest = TRUE, dig.lab = 4)

# ---- Figures --------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))
mp <- function(col, title, leg, dat = G)
  ggplot(dat, aes(X_UTM33N, Y_UTM33N, colour = .data[[col]])) +
  geom_point(size = 0.45) + coord_equal() +
  labs(title = title, colour = leg, x = NULL, y = NULL) +
  theme(axis.text = element_blank())

# F20: operational map of the thickness expected by the forest
g1 <- mp("class_rf", "Expected thickness (random forest, ensemble)", "class (cm)") +
      scale_colour_viridis_d(option = "turbo")
ggplot2::ggsave(file.path(PATHS$figs, "F20_map_rf.png"), g1,
                width = 7, height = 6, dpi = 200)

# F23: e2tree interpretation = rule zones (NOT a prediction)
g4 <- mp("rule_zone",
         "e2tree interpretation: rule zones of the ensemble", "node/rule") +
      scale_colour_viridis_d(option = "turbo")
ggplot2::ggsave(file.path(PATHS$figs, "F23_map_rules.png"), g4,
                width = 7, height = 6, dpi = 200)

# F24: reliability (training domain vs extrapolation)
g5 <- ggplot(G, aes(X_UTM33N, Y_UTM33N,
                     colour = ifelse(in_domain, "in domain", "extrapolation"))) +
  geom_point(size = 0.45) + coord_equal() +
  scale_colour_manual(values = c("extrapolation" = "#d73027",
                                 "in domain" = "#4575b4")) +
  labs(title = "Reliability: training domain vs extrapolation",
       colour = NULL, x = NULL, y = NULL) +
  theme(axis.text = element_blank())
ggplot2::ggsave(file.path(PATHS$figs, "F24_map_extrapolation.png"), g5,
                width = 7, height = 6, dpi = 200)

# F25: operational map = ensemble only where reliable (in domain)
Gin <- G[G$in_domain, ]
g6 <- mp("class_rf",
         "Expected thickness (ensemble) - reliable areas only (in domain)",
         "class (cm)", dat = Gin) + scale_colour_viridis_d(option = "turbo")
ggplot2::ggsave(file.path(PATHS$figs, "F25_map_rf_indomain.png"), g6,
                width = 7, height = 6, dpi = 200)

out <- G[, c("X_UTM33N", "Y_UTM33N", "pred_rf", "class_rf",
             "rule_zone", "in_domain", "extrap_vars")]
names(out)[names(out) == "pred_rf"] <- "pred_ensemble_RF"
write.csv(out, file.path(PATHS$tables, "T12_grid_predictions.csv"),
          row.names = FALSE)
saveRDS(G, file.path(PATHS$models, "grid_pred.rds"))
message("05_predict_maps OK | prediction=ensemble RF | e2tree=interpretation | ",
        "in-domain=", sum(G$in_domain), "/", nrow(G))
