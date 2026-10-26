# 13_figures_en.R - Manuscript figures
# Rebuilt from the saved results, without refitting any model: run it after
# the rest of the pipeline.
source("R/00_setup.R")
suppressPackageStartupMessages({library(ggplot2); library(e2tree); library(ranger); library(rpart.plot)})
TB <- PATHS$tables; FG <- PATHS$figs
rd <- function(f) read.csv(file.path(TB, f), check.names = FALSE)
theme_set(theme_minimal(base_size = 11))
# Charts as vector PDF and 300 dpi PNG, maps as 600 dpi PNG (save_gg,
# save_base and save_map in 00_setup.R).
sv <- function(n, p, w = 7, h = 5, vector = TRUE)
  save_gg(n, p, w, h, vector = vector)

# ---- F30: accuracy and structural fidelity (nLoI) --------------------------
# The values come from T15_accuracy_fidelity_tradeoff.csv, computed by
# R/03c_tradeoff_nloi.R, so figure and table cannot diverge.
to <- rd("T15_accuracy_fidelity_tradeoff.csv")
to$config <- sub("RF predictive \\(nt1200,ns5\\)", "accuracy-optimal (nt1200, ns5)", to$config)
to$config <- sub("nt600 ns25", "nt600, ns25", to$config)
to$config <- sub("nt400 ns50", "nt400, ns50", to$config)
to$config <- sub("nt300 ns80 \\[selected\\]", "nt300, ns80 (selected)", to$config)
to$isP <- grepl("accuracy-optimal", to$config); to$isC <- grepl("selected", to$config)
# Red and green from COL_REGIME (00_setup.R), distinguishable also by readers
# who cannot tell the two apart; shapes and labels give a second cue.
col30 <- ifelse(to$isP, COL_REGIME[["red"]],
                ifelse(to$isC, COL_REGIME[["green"]], "#2c7fb8"))
save_base("F30_tradeoff", function() {
  op <- par(mar = c(5, 5.5, 4, 2), xpd = NA)
  on.exit(par(op))
  plot(to$RF_RMSE, to$nLoI, type = "n",
       xlim = c(min(to$RF_RMSE) - 0.6, max(to$RF_RMSE) + 0.9),
       ylim = c(max(to$nLoI) * 1.15, 0),
       xlab = "RMSE (cm)",
       ylab = "nLoI",
       main = "")
  lines(to$RF_RMSE[order(to$RF_RMSE)], to$nLoI[order(to$RF_RMSE)],
        col = "grey60", lty = 2)
  points(to$RF_RMSE, to$nLoI,
         pch = ifelse(to$isP, 4, ifelse(to$isC, 17, 19)),
         col = col30, cex = ifelse(to$isC, 2.4, 1.8), lwd = 2)
  text(to$RF_RMSE, to$nLoI, to$config, pos = 3, cex = 0.72, col = col30)
  legend("bottomright", bty = "n", cex = 0.8, pch = c(4, 19, 17),
         col = c(COL_REGIME[["red"]], "#2c7fb8", COL_REGIME[["green"]]),
         legend = c("accuracy-optimal RF", "regularised RF",
                    "selected RF"))
}, w = 1200 / 150, h = 800 / 150)

# ---- F40: curve of the local data required -----------------------------------
# In the foreground, as in the text, the curve of the delivered model
# (covariates only: solid line and points); RFsp and RFsp+RK are the
# interpolation references (dashed, thin). The operational threshold
# R2 = 0.40 is marked: the delivered curve crosses it at 25% of in-area data.
cu <- rd("T23_data_requirement_curve.csv")
lv40 <- c("RF (delivered)", "RFsp (benchmark)",
          "RFsp + RK (benchmark)")

# Secondary axis in absolute units (R/31_curve_density.R). The density is
# linear in p up to the per-block floor() (difference < 1%), so the conversion
# factor is taken from the end of the curve.
sec <- NULL
if (file.exists(file.path(TB, "T44_curve_density.csv"))) {
  den <- rd("T44_curve_density.csv")
  kdn <- den$dens_pooled[which.max(den$p)] / (max(den$p) * 100)
  # The ticks of the upper axis fall on the same percentages as the lower one
  # (0, 10, ..., 50%): each density is read above its percentage.
  sec <- sec_axis(~ . * kdn,
                  name = expression("In-area data (pt/km"^2*", pooled)"),
                  breaks = seq(0, 50, 10) * kdn,
                  labels = function(b) formatC(b, format = "f", digits = 2))
}

# The two thresholds and the bands of the three regimes, with the colours of
# Figure 4 (COL_REGIME) strongly lightened, so the curve is read against the
# regime it falls in.
y40 <- range(c(cu$R2_rf, cu$R2_rfsp, cu$R2_rfsp_rk, 0.20, 0.40))
y40 <- c(min(0.10, y40[1] - 0.03), y40[2] + 0.03)
p40 <- ggplot(cu, aes(p * 100)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.20,
           fill = COL_REGIME[["red"]], alpha = 0.10) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.20, ymax = 0.40,
           fill = COL_REGIME[["amber"]], alpha = 0.14) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.40, ymax = Inf,
           fill = COL_REGIME[["green"]], alpha = 0.10) +
  geom_hline(yintercept = c(0.20, 0.40), linetype = 3, colour = "grey35") +
  # labels on the right, just above the thresholds (the curves are above 0.45
  # there); thresholds as strings, because plotmath would print 0.40 as "0.4"
  annotate("text", x = 50, y = c(0.405, 0.205, y40[1] + 0.005), hjust = 1,
           vjust = 0, size = 3, colour = "grey30", parse = TRUE,
           label = c("'operational use ('*R^2 >= '0.40'*')'",
                     "'in-area support ('*'0.20' <= R^2*' < 0.40)'",
                     "'not usable ('*R^2 < '0.20'*')'")) +
  geom_line(aes(y = R2_rfsp_rk, colour = "RFsp + RK (benchmark)"),
            linewidth = 0.6, linetype = 3) +
  geom_line(aes(y = R2_rfsp, colour = "RFsp (benchmark)"),
            linewidth = 0.5, linetype = 2) +
  geom_line(aes(y = R2_rf, colour = "RF (delivered)"),
            linewidth = 1.2) +
  geom_point(aes(y = R2_rf, colour = "RF (delivered)"),
             size = 2.8) +
  scale_colour_manual(values = c("RF (delivered)" = "#e6550d",
                                 "RFsp (benchmark)" = "#2c7fb8",
                                 "RFsp + RK (benchmark)" = "#1a9850"),
                      breaks = lv40) +
  scale_x_continuous(breaks = seq(0, 50, 10),
                     sec.axis = if (is.null(sec)) waiver() else sec) +
  coord_cartesian(ylim = y40) +
  labs(x = "In-area data (%)",
       y = expression(R^2~"(spatial cross-validation)"), colour = NULL) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = NA))
sv("F40_data_requirement.png", p40, 8, 5.2)

# ---- F41: operating regimes, from the scenario registry ----------------------
# Source: T47_scenario_registry.csv (R/30_scenario_registry.R), with the 11
# designs of Table 4, the same acronyms and the same order. Thresholds as in
# the text (green >= 0.40, amber 0.20-0.40, red < 0.20); R5 apart, as the
# reference with leakage. Crosses: raw GPT on the same folds, where it exists.
va <- rd("T47_scenario_registry.csv")
va$cls <- ifelse(va$acr == "R5", "lk",
          ifelse(va$RFsp >= 0.40, "op", ifelse(va$RFsp >= 0.20, "ia", "nu")))
va$GPT <- ifelse(is.na(va$GPT_R2), va$R2_GPT, va$GPT_R2)
va$acr <- factor(va$acr, levels = rev(va$acr))   # first design at the top
p41 <- ggplot(va, aes(acr, RFsp, fill = cls)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(0.20, 0.40), linetype = "dashed",
             colour = "grey40", linewidth = 0.35) +
  geom_point(aes(y = GPT), shape = 4, size = 3, colour = "black",
             show.legend = FALSE, na.rm = TRUE) +
  # COL_REGIME (00_setup.R): green, amber and red as the text calls them,
  # distinguishable also by readers who cannot tell red from green.
  scale_fill_manual(
    values = c(op = COL_REGIME[["green"]], ia = COL_REGIME[["amber"]],
               nu = COL_REGIME[["red"]], lk = "grey65"),
    breaks = c("op", "ia", "nu", "lk"),
    # thresholds as strings: plotmath would print 0.40 as "0.4"
    labels = c(expression(R^2 >= "0.40"~"(operational use)"),
               expression("0.20" <= R^2~"< 0.40 (in-area support)"),
               expression(R^2 < "0.20"~"(not usable)"),
               expression("R5: random CV (leakage reference)")),
    name = NULL) +
  labs(x = "Validation design", y = expression("Cross-validated"~R^2~"(RFsp)")) +
  theme(legend.position = "bottom",
        legend.text = element_text(hjust = 0)) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
sv("F41_validity_domain.png", p41, 8, 4.8)

# ---- F46: comparison with the empirical models -------------------------------
# Acronyms and order from the registry; the x axis does not quote table
# numbers, which change with the layout.
be  <- rd("T28_fair_benchmark_full.csv")
reg <- rd("T47_scenario_registry.csv")
lab <- setNames(reg$acr, reg$key)
L <- data.frame(
  scenario = factor(lab[be$scenario], levels = intersect(reg$acr, lab[be$scenario])),
  RFsp = be$RFsp_R2, GPT = be$GPT_R2, SAPT = be$SAPT_R2, SEPT = be$SEPT_R2)
Lg <- reshape(L, varying = c("RFsp","GPT","SAPT","SEPT"), v.names = "R2",
              timevar = "model", times = c("RFsp","GPT","SAPT","SEPT"),
              direction = "long")
# legend in the order of the columns of Table 8: RFsp, then the empirical models
Lg$model <- factor(Lg$model, levels = c("RFsp", "GPT", "SAPT", "SEPT"))
# RFsp in blue, as in Figure 5; the three empirical models in greys from
# darkest to lightest, with a dark outline that keeps the light grey visible
# on the white background. No red or green, which are reserved for the
# regimes of Figure 4.
p46 <- ggplot(Lg, aes(scenario, R2, fill = model)) +
  geom_col(position = "dodge", colour = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c(RFsp = "#2c7fb8", GPT = "#404040",
                               SAPT = "#A6A6A6", SEPT = "#D9D9D9"),
                    labels = c(GPT="GPT (raw)", RFsp="RFsp",
                               SAPT="SAPT (raw)", SEPT="SEPT (raw)")) +
  labs(x = "Validation design", y = expression("Cross-validated"~R^2),
       fill = NULL)
sv("F46_fair_benchmark_full.png", p46, 9, 5)

# ---- F11: native measures of structural fidelity ----------------------------
me <- rd("T6_e2tree_fidelity.csv")
# Only the measures defined in the text (no wRMSE); the permutation reference
# quoted in the caption is drawn with crosses.
m <- me[!grepl("wRMSE", me$method, ignore.case = TRUE), ]
m <- m[order(m$observed), ]
pm <- if ("perm_mean" %in% names(m)) m$perm_mean else
      if ("null_mean" %in% names(m)) m$null_mean else NULL
# Vector PDF and PNG (save_base).
save_base("F11_fidelity_measures", function() {
  op <- par(mar = c(8, 4.5, 2, 1))
  on.exit(par(op))
  bp <- barplot(m$observed, names.arg = m$method, las = 2, col = "#2c7fb8",
                ylab = "observed value", ylim = c(0, max(
                  c(m$observed, m$ci_upper, m$perm_mean), na.rm = TRUE) * 1.15))
  if (!is.null(m$ci_lower))
    arrows(bp, m$ci_lower, bp, m$ci_upper, angle = 90, code = 3,
           length = 0.04, col = "grey30")
  if (!is.null(pm)) {
    points(bp, pm, pch = 4, cex = 1.3, lwd = 2, col = "#d73027")
    legend("topleft", bty = "n", cex = 0.85,
           legend = c("observed (95% CI)", "permutation null (mean)"),
           pch = c(15, 4), col = c("#2c7fb8", "#d73027"))
  }
}, w = 1150 / 150, h = 720 / 150)

# ---- F10: explanatory tree ---------------------------------------------------
e2 <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))
# Node labels renamed for the figure only: ThickIso becomes z0 (the initial
# thickness, not z), Elevation becomes Altitude. The column names do not
# change: only the copy passed to the plot is modified.
e2_disp <- function(obj) {
  ren <- c(ThickIso = "z0", Elevation = "Altitude")
  tr <- obj$tree
  for (cl in c("variable", "splitLabel", "variableSur", "splitLabelSur", "path"))
    if (cl %in% names(tr$tree)) {
      v <- as.character(tr$tree[[cl]])
      for (o in names(ren)) v <- gsub(o, ren[[o]], v, fixed = TRUE)
      tr$tree[[cl]] <- v
    }
  if (!is.null(tr$data))
    names(tr$data)[match(names(ren), names(tr$data), nomatch = 0)] <-
      unname(ren)[names(ren) %in% names(tr$data)]
  tr
}
# plot_e2tree() fixes type = 4, extra = 101, cex = 0.9, box.palette = "Blues"
# and shadow.col = "gray" in its own signature, so they cannot be changed
# through `...` ("matched by multiple actual arguments"); with 23 leaves those
# values make the labels overlap. rpart.plot is therefore used on the
# package's converter:
#   type = 2, extra = 0: the figure shows the structure and the node numbers;
#     sample sizes and grid shares are in Table S5;
#   nn = TRUE: the text cites node 105 and Figure 9 uses the same numbers;
#   width equal to the \textwidth of the template (468.33 pt, i.e. 6.48 in):
#     the figure enters at scale 1:1 and the text has the intended size.
TEXTWIDTH_IN <- 468.33 / 72.27
f10_plot <- function(tr) {
  save_base("F10_e2tree", function() {
    par(mar = c(0.1, 0.1, 0.1, 0.1))  # par settings die with the device: no restore needed
    rpart.plot::rpart.plot(
      e2tree::rpart2Tree(tr, e2$rf), main = "", type = 2, extra = 0,
      nn = TRUE, box.palette = "Grays", shadow.col = NA, cex = NULL,
      tweak = 1.10, faclen = 0, varlen = 0, roundint = FALSE)
  }, w = TEXTWIDTH_IN, h = 3.5)
}
tryCatch(f10_plot(e2_disp(e2)), error = function(e) {
  message("F10 (renamed) failed: ", conditionMessage(e),
          " -> falling back to the original names")
  f10_plot(e2$tree)
})

# ---- F42-F45: operational maps, from the saved grid ---------------------------
# Shared functions of 00_setup.R (gg_field_map, gg_class_map): geom_raster and
# perceptually uniform palettes, volcanoes, scale bar and north arrow. The
# isopachs are computed on a field smoothed for display only; the q95 classes
# have a modal filter that never lowers the class of a cell.
G <- rd("T27_operational_grid.csv")
G$in_domain <- as.logical(G$in_domain)
# The operational map comes from a forest with the covariates only: it has no
# spatial-extrapolation artefacts and is shown over the whole grid, without a
# mask. Volcano labels at text_size 3.3, with Roccamonfina; the panel of the
# interval width carries the values on the contours.
#  - every panel has a title that says which quantity it shows;
#  - the outline of the sampled zone (cells within the variogram range,
#    1,206 m, of a measurement point; 38_applicability.R): cyan with a black
#    edge, so it is not confused with the contours (white and black) and
#    stays visible on every palette. Outside the outline no training point
#    falls within the range, as in B15b and NNDM: that is the red regime.
zone <- readRDS(file.path(PATHS$models, "sampled_zone.rds"))$zone
zdf  <- as.data.frame(sf::st_coordinates(zone))
zdf$grp <- interaction(zdf$L1, zdf$L2, if ("L3" %in% names(zdf)) zdf$L3 else 1)
zone_layers <- list(
  geom_path(data = zdf, aes(X, Y, group = grp), inherit.aes = FALSE,
            colour = "black", linewidth = 1.1),
  geom_path(data = zdf, aes(X, Y, group = grp), inherit.aes = FALSE,
            colour = "#00E5FF", linewidth = 0.55))
map_title <- theme(plot.title = element_text(size = 11, face = "bold"))
save_map("F42_operational_map.png",
  gg_field_map(G, "pred", "Expected thickness", "cm", palette = "viridis",
               isopach = ISOPACH_CM, carto = TRUE,
               carto_args = list(text_size = 3.3, roccamonfina = TRUE,
                                 cilento = TRUE)) + zone_layers + map_title)
save_map("F43_uncertainty.png",
  gg_field_map(G, "pi_width_smooth", "Width of the 90% prediction interval",
               "cm", palette = "inferno", direction = -1,
               isopach = UNC_CONTOURS, carto = TRUE,
               carto_args = list(text_size = 3.3, roccamonfina = TRUE,
                                 cilento = TRUE)) + zone_layers + map_title)
# Reversed legend: the highest class at the top.
save_map("F44_conservative_q95.png",
  gg_class_map(G, "q95", Q95_BREAKS, Q95_LABELS,
               "0.95 conditional quantile of thickness", "cm",
               modal_w = 3, carto = TRUE, reverse_legend = TRUE,
               carto_args = list(text_size = 3.3, roccamonfina = TRUE,
                                 cilento = TRUE)) + zone_layers + map_title)
save_map("F45_rule_zones.png",
  ggplot(G, aes(X_UTM33N, Y_UTM33N, fill = factor(rule_zone))) +
    geom_raster() + coord_equal() +
    scale_fill_viridis_d(option = "turbo",
                         guide = guide_legend(ncol = 2)) +
    labs(title = "Rule zones of the explanatory tree", fill = "rule",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    carto_layers(G, text_size = 3.3, roccamonfina = TRUE, cilento = TRUE) +
    zone_layers + map_title)

# ---- F58: RFsp and the delivered forest, on the same colour scale -----------
# The figure of section fullcov: same raster, two models. In the Cilento
# (southern band, no training points) RFsp produces a spurious cover, whereas
# the covariates-only forest gives the expected thin values. With facet_wrap
# the panels are identical and there is a single legend; an ellipse and a
# label mark the Cilento.
dcl <- readRDS(file.path(PATHS$models, "data_clean.rds"))
opm <- readRDS(file.path(PATHS$models, "operational.rds"))
# Covariates from predictor_vars() (00_setup.R): the 16 variables of
# Table 2.
pv58 <- predictor_vars(dcl$model_df)
g58 <- dcl$newdata[stats::complete.cases(
  dcl$newdata[, c(pv58, "X_UTM33N", "Y_UTM33N")]), ]
g58$rfsp <- predict(opm$rfsp, g58)$predictions
M58 <- merge(g58[, c("X_UTM33N", "Y_UTM33N", "rfsp")],
             G[, c("X_UTM33N", "Y_UTM33N", "pred")],
             by = c("X_UTM33N", "Y_UTM33N"))
DD <- rbind(
  data.frame(X_UTM33N = M58$X_UTM33N, Y_UTM33N = M58$Y_UTM33N, v = M58$rfsp,
             panel = "(a)"),
  data.frame(X_UTM33N = M58$X_UTM33N, Y_UTM33N = M58$Y_UTM33N, v = M58$pred,
             panel = "(b)"))
th58 <- seq(0, 2 * pi, length.out = 240)
ell58 <- data.frame(x = CILENTO_ELLIPSE$cx + CILENTO_ELLIPSE$ax * cos(th58),
                    y = CILENTO_ELLIPSE$cy + CILENTO_ELLIPSE$ay * sin(th58))
p58 <- ggplot(DD, aes(X_UTM33N, Y_UTM33N, fill = v)) +
  geom_raster() + coord_equal() +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_viridis_c(option = "viridis") +
  geom_path(data = ell58, aes(x, y), inherit.aes = FALSE,
            colour = "white", linewidth = 0.7) +
  annotate("text", x = CILENTO_ELLIPSE$cx, y = CILENTO_ELLIPSE$cy + 26500,
           label = "Cilento",
           colour = "white", size = 3.4, fontface = "bold") +
  labs(fill = "cm", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(hjust = 0, size = 11)) +
  carto_layers(DD, north = TRUE, text_size = 3.0, roccamonfina = TRUE)
sv("F58_covonly_vs_rfsp.png", p58, 11, 5, vector = FALSE)   # map: 600 dpi PNG

message("13_figures_en OK | manuscript figures regenerated")
