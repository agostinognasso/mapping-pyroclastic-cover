# 00_setup.R - Global settings: paths, seed, libraries, shared helpers.
# Part of the analysis for "Spatial cross-validation and sampling-requirement
# curves for mapping pyroclastic cover thickness" (random forests, spatial
# cross-validation and an explainable ensemble tree, Campania, Italy).

# ---- Paths ----------------------------------------------------------------
# Run every script from the repository root (the folder that contains data/
# and R/). Input data live in data/.
PROJ <- getwd()
DATA <- file.path(PROJ, "data")
if (!file.exists(file.path(DATA, "dataset.csv")))
  stop("data/dataset.csv not found. Run the scripts from the repository root ",
       "(the folder that contains data/ and R/). Current working directory: ", PROJ)

# Output folders. By default models, CSV tables and GeoTIFFs go to output/,
# figures to output/figs/ and the LaTeX tables to output/tex/. To run the
# pipeline without overwriting existing results, set the environment variables
# below, for example:
#   PTCAM_OUT=output_new PTCAM_PAPER=output_new/tex PTCAM_FIGS=output_new/figs Rscript run_all.R
OUTDIR   <- file.path(PROJ, Sys.getenv("PTCAM_OUT",   unset = "output"))
PAPERDIR <- file.path(PROJ, Sys.getenv("PTCAM_PAPER", unset = "output/tex"))
FIGDIR   <- file.path(PROJ, Sys.getenv("PTCAM_FIGS",  unset = "output/figs"))

PATHS <- list(
  root    = PROJ,
  data    = file.path(DATA, "dataset.csv"),
  newdata = file.path(DATA, "prediction_grid.csv"),
  out     = OUTDIR,
  paper   = PAPERDIR,
  figs    = FIGDIR,
  tables  = file.path(OUTDIR, "tables"),
  models  = file.path(OUTDIR, "models")
)
for (d in c(PATHS$figs, PATHS$tables, PATHS$models, PATHS$paper))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Global seed ----------------------------------------------------------
SEED <- 12345
set.seed(SEED)

# ---- Libraries ------------------------------------------------------------
suppressPackageStartupMessages({
  library(ranger)
  library(randomForest)
  library(e2tree)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ---- Variable definitions -------------------------------------------------
TARGET <- "thick_cm"

# Identifier columns ("localitiy" is the column name as published in the
# PT-Cam-2 data).
ID_COLS        <- c("X_UTM33N", "Y_UTM33N", "note", "localitiy", "type")
# Empirical thickness models used only as benchmarks, never as covariates
EMPIRICAL_COLS <- c("GPT_model", "SEPT_model", "SAPT_model")
# Columns excluded from the covariates. DistSource would give unrealistic
# predictions. DenHydnet is not among the covariates documented in Ebrahimi
# et al. (2024): its origin is unknown, it is zero in 7167 of 8109 records
# (88.4%), and it is the weakest covariate (%IncMSE 6.3 against 20.0 for the
# next weakest). Without the two columns there are 16 covariates, those of
# Table 2 of the paper.
EXCL_RF        <- c("DistSource", "DenHydnet")

# Covariates: every column except identifiers, response, empirical models and
# the excluded columns
predictor_vars <- function(df) {
  setdiff(names(df), c(TARGET, ID_COLS, EMPIRICAL_COLS, EXCL_RF))
}

# ---- Metrics --------------------------------------------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))
mae  <- function(obs, pred) mean(abs(obs - pred))
r2   <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
bias <- function(obs, pred) mean(pred - obs)

metrics_tbl <- function(obs, pred, model = "model") {
  data.frame(model = model,
             RMSE = rmse(obs, pred), MAE = mae(obs, pred),
             R2 = r2(obs, pred), Bias = bias(obs, pred))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Assign each observation to a terminal node of the e2tree ---------------
# Evaluates the path of each rule (a conjunction of conditions) on `data`.
# The paths partition the space: each observation goes to the first rule it
# satisfies, which is deterministic even with surrogate splits.
e2_assign <- function(tree, data) {
  tn    <- e2tree::nodes(tree, terminal = TRUE)
  paths <- as.character(tn$path)
  node  <- as.integer(as.character(tn$node))
  predv <- as.numeric(as.character(tn$pred))
  n <- nrow(data)
  out_node <- rep(NA_integer_, n); out_pred <- rep(NA_real_, n)
  to_expr <- function(p) {
    parts <- trimws(strsplit(p, " & ", fixed = TRUE)[[1]])
    parts <- ifelse(startsWith(parts, "!"),
                     paste0("!(", sub("^!\\s*", "", parts), ")"),
                     paste0("(", parts, ")"))
    paste(parts, collapse = " & ")
  }
  for (k in seq_along(paths)) {
    if (is.na(paths[k]) || !nzchar(paths[k])) {
      mask <- rep(TRUE, n)
    } else {
      ex <- tryCatch(eval(parse(text = to_expr(paths[k])),
                          envir = as.list(data)),
                      error = function(e) rep(FALSE, n))
      mask <- !is.na(ex) & ex
    }
    sel <- mask & is.na(out_node)
    out_node[sel] <- node[k]; out_pred[sel] <- predv[k]
  }
  data.frame(node = out_node, pred_node = out_pred)
}

# ---- Maps: shared rendering functions ---------------------------------------
# Every continuous-field map (thickness, q95, uncertainty) goes through here:
#  - geom_raster on the real cell, not geom_point on the grid;
#  - perceptually uniform palettes (viridis, mako, inferno): turbo has
#    luminance reversals that draw bands that are not in the data;
#  - isopach: given a vector of levels in cm, overlays labelled isopachs
#    computed with terra::as.contour.
# A single definition keeps all scripts drawing maps the same way.

# Isopachs from a regular grid: lines and labels for ggplot.
#  - smooth_w: side, in cells, of a moving average applied only to the field
#    the contours are drawn from; the fill stays the raw field. Cell-by-cell
#    predictions give fragmented contours, and the average makes them
#    readable without touching the data;
#  - min_pts: drops contour pieces with fewer vertices (fragments and noise
#    rings), so labels fall on real lines.
field_contours <- function(dd, col, levels, crs = "EPSG:32633",
                           smooth_w = NULL, min_pts = 0) {
  if (!requireNamespace("terra", quietly = TRUE))
    stop("field_contours: the 'terra' package is required")
  r  <- terra::rast(data.frame(x = dd$X_UTM33N, y = dd$Y_UTM33N, z = dd[[col]]),
                    type = "xyz", crs = crs)
  if (!is.null(smooth_w))
    r <- terra::mask(terra::focal(r, w = smooth_w, fun = "mean",
                                  na.rm = TRUE, na.policy = "omit"), r)
  cn <- terra::as.contour(r, levels = levels)
  lv <- terra::values(cn)$level                # level of each geometry
  g  <- terra::geom(cn)                         # geom, part, x, y, hole
  lines <- data.frame(x = g[, "x"], y = g[, "y"],
                      grp = paste(g[, "geom"], g[, "part"], sep = "_"),
                      level = lv[g[, "geom"]])
  if (min_pts > 0) {
    keep  <- names(which(table(lines$grp) >= min_pts))
    lines <- lines[lines$grp %in% keep, ]
  }
  if (nrow(lines) == 0)
    return(list(lines = lines,
                labels = data.frame(x = numeric(0), y = numeric(0),
                                    level = numeric(0))))
  # One label per level, at the right-most vertex of the longest piece: with
  # nearly concentric isopachs the labels separate along x, whereas on the
  # centroids they would overlap.
  labs <- do.call(rbind, lapply(split(lines, lines$level), function(s) {
    parts   <- split(s, s$grp)
    longest <- parts[[which.max(vapply(parts, nrow, integer(1)))]]
    i <- which.max(longest$x)
    data.frame(x = longest$x[i], y = longest$y[i], level = longest$level[i])
  }))
  list(lines = lines, labels = labs)
}

# ---- Shared cartographic elements (volcanoes, scale bar, north arrow) -------
# UTM 33N (EPSG:32633) coordinates projected from WGS84:
#   Somma-Vesuvio  40.8214 N, 14.4260 E  ->  (451596, 4519090)
#   Campi Flegrei  40.8275 N, 14.1391 E  ->  (427408, 4519965)
VENTS_UTM33 <- data.frame(
  name = c("Vesuvius", "Campi Flegrei"),
  X = c(451596, 427408), Y = c(4519090, 4519965))

# Roccamonfina (41.290 N, 13.977 E -> 414342, 4571455): the third volcanic
# district in the study area. Its fallout contributes little to z0 in the
# north of the area, so it is drawn like the other volcanoes (white triangle,
# label on a white background). It is a separate data.frame because maps show
# it only with roccamonfina = TRUE.
ROCCAMONFINA_UTM33 <- data.frame(name = "Roccamonfina",
                                 X = 414342, Y = 4571455)

# Southern Cilento: a distal area without measurements, used in the paper to
# check the delivered map. The same ellipse bounds the area in Figure S10
# (13_figures_en.R) and in the numbers quoted in the text (34_cilento_check.R):
# centre (523500, 4474000) m, semi-axes 28 and 21 km.
CILENTO_ELLIPSE <- list(cx = 523500, cy = 4474000, ax = 28000, ay = 21000)
# label 6 km north of the centre: at the centre it touches the 50 cm isopach
CILENTO_UTM33   <- data.frame(name = "Cilento", X = CILENTO_ELLIPSE$cx,
                              Y = CILENTO_ELLIPSE$cy + 6000)
in_cilento <- function(x, y)
  ((x - CILENTO_ELLIPSE$cx) / CILENTO_ELLIPSE$ax)^2 +
  ((y - CILENTO_ELLIPSE$cy) / CILENTO_ELLIPSE$ay)^2 <= 1

# Reusable ggplot layers: volcano symbols, scale bar and north arrow, in the
# south-west corner of the grid, which falls in the Tyrrhenian Sea and covers
# no cells. The grid is in UTM, with north up.
carto_layers <- function(dd, vents = TRUE, scalebar_km = 20, north = TRUE,
                         text_size = 2.7, roccamonfina = FALSE,
                         cilento = FALSE) {
  xr <- range(dd$X_UTM33N, na.rm = TRUE); yr <- range(dd$Y_UTM33N, na.rm = TRUE)
  x0 <- xr[1] + 2000; y0 <- yr[1] + 5000
  ly <- list()
  if (isTRUE(roccamonfina)) ly <- c(ly, list(
    ggplot2::geom_point(data = ROCCAMONFINA_UTM33, ggplot2::aes(X, Y),
                        inherit.aes = FALSE, shape = 24, fill = "white",
                        colour = "black", size = 2.4, stroke = 0.6),
    ggplot2::geom_label(data = ROCCAMONFINA_UTM33, ggplot2::aes(X, Y, label = name),
                        inherit.aes = FALSE, size = text_size,
                        colour = "grey10", fill = "white", alpha = 0.75,
                        linewidth = 0, vjust = -0.55)))
  # Name of an area, not of a point: italic and without a symbol
  # (cartographic convention), centred on the Cilento ellipse.
  if (isTRUE(cilento)) ly <- c(ly, list(
    ggplot2::geom_label(data = CILENTO_UTM33, ggplot2::aes(X, Y, label = name),
                        inherit.aes = FALSE, size = text_size,
                        colour = "grey10", fill = "white", alpha = 0.75,
                        linewidth = 0, fontface = "italic")))
  if (isTRUE(vents)) ly <- c(ly, list(
    ggplot2::geom_point(data = VENTS_UTM33, ggplot2::aes(X, Y),
                        inherit.aes = FALSE, shape = 24, fill = "white",
                        colour = "black", size = 2.4, stroke = 0.6),
    ggplot2::geom_label(data = VENTS_UTM33, ggplot2::aes(X, Y, label = name),
                        inherit.aes = FALSE, size = text_size,
                        colour = "grey10", fill = "white", alpha = 0.75,
                        linewidth = 0, vjust = -0.55)))
  if (!is.null(scalebar_km)) {
    L <- scalebar_km * 1000
    ly <- c(ly, list(
      ggplot2::annotate("segment", x = x0, xend = x0 + L, y = y0, yend = y0,
                        linewidth = 0.9, colour = "grey15"),
      ggplot2::annotate("segment", x = c(x0, x0 + L), xend = c(x0, x0 + L),
                        y = y0 - 1500, yend = y0 + 1500,
                        linewidth = 0.9, colour = "grey15"),
      ggplot2::annotate("text", x = x0 + L / 2, y = y0 + 4200,
                        label = paste(scalebar_km, "km"),
                        size = text_size, colour = "grey15")))
  }
  if (isTRUE(north)) {
    xn <- x0 + 4000; yn <- y0 + 12000
    ly <- c(ly, list(
      ggplot2::annotate("segment", x = xn, xend = xn, y = yn, yend = yn + 9000,
                        linewidth = 0.7, colour = "grey15",
                        arrow = grid::arrow(length = grid::unit(5, "pt"),
                                            type = "closed")),
      ggplot2::annotate("text", x = xn, y = yn + 13000, label = "N",
                        size = text_size + 0.5, fontface = "bold",
                        colour = "grey15")))
  }
  ly
}

# Continuous-field map with optional isopachs; all maps of the paper use it.
#  - title = NULL: no title in the figure, the caption carries it;
#  - carto = TRUE: volcanoes, scale bar and north arrow (carto_layers);
#  - isopach_smooth, contour_min_pts: smoothed isopachs without fragments
#    (see field_contours); the fill stays the raw field.
gg_field_map <- function(dd, col, title = NULL, legend = "cm",
                         palette = "viridis", direction = 1,
                         limits = NULL, isopach = NULL, label_size = 3.1,
                         contour_light = NULL, carto = FALSE,
                         carto_args = list(),
                         isopach_smooth = 5, contour_min_pts = 10) {
  g <- ggplot2::ggplot(dd, ggplot2::aes(X_UTM33N, Y_UTM33N)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data[[col]])) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_viridis_c(option = palette, direction = direction,
                                  limits = limits) +
    ggplot2::labs(title = title, fill = legend, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  if (!is.null(isopach)) {
    ct <- field_contours(dd, col, isopach,
                         smooth_w = isopach_smooth, min_pts = contour_min_pts)
    # Double stroke, white below and dark above: the isopachs stay readable on
    # both the dark background and the bright end of viridis.
    g <- g +
      ggplot2::geom_path(data = ct$lines, ggplot2::aes(x, y, group = grp),
                         colour = "white", linewidth = 0.6, alpha = 0.7,
                         inherit.aes = FALSE) +
      ggplot2::geom_path(data = ct$lines, ggplot2::aes(x, y, group = grp),
                         colour = "grey10", linewidth = 0.25,
                         inherit.aes = FALSE) +
      ggplot2::geom_label(data = ct$labels, ggplot2::aes(x, y, label = level),
                          inherit.aes = FALSE, size = label_size,
                          colour = "grey10", fill = "white", alpha = 0.8,
                          linewidth = 0)
  }
  # Light contours without labels (for example in F43): thin, semi-transparent
  # white lines that give a visual reference without suggesting a precision
  # that an already smoothed field does not have.
  if (!is.null(contour_light)) {
    ctl <- field_contours(dd, col, contour_light, min_pts = contour_min_pts)
    g <- g + ggplot2::geom_path(data = ctl$lines, ggplot2::aes(x, y, group = grp),
                                colour = "white", linewidth = 0.25, alpha = 0.5,
                                inherit.aes = FALSE)
  }
  if (isTRUE(carto)) g <- g + do.call(carto_layers, c(list(dd), carto_args))
  g
}

# Standard isopach levels for the expected thickness (cm).
ISOPACH_CM <- c(25, 50, 100, 150, 200, 300)

# Map of discrete classes, for the conservative scenario read as a hazard
# zoning. modal_w (for example 3) applies a modal filter to the classes: the
# final class is the maximum of the original one and the local mode, so
# isolated cells are corrected only upwards and no cell is downgraded.
gg_class_map <- function(dd, col, breaks, labels, title = NULL, legend = "cm",
                         palette = "inferno", direction = -1,
                         modal_w = NULL, carto = FALSE,
                         carto_args = list(), reverse_legend = FALSE) {
  dd2 <- dd
  ci  <- as.integer(cut(dd[[col]], breaks = breaks, labels = FALSE,
                        include.lowest = TRUE, right = FALSE))
  if (!is.null(modal_w)) {
    if (!requireNamespace("terra", quietly = TRUE))
      stop("gg_class_map: the 'terra' package is required for modal_w")
    r  <- terra::rast(data.frame(x = dd$X_UTM33N, y = dd$Y_UTM33N, z = ci),
                      type = "xyz", crs = "EPSG:32633")
    rm_ <- terra::mask(terra::focal(r, w = modal_w, fun = "modal",
                                    na.rm = TRUE, na.policy = "omit"), r)
    sm <- terra::extract(rm_, cbind(dd$X_UTM33N, dd$Y_UTM33N))[, 1]
    ci <- pmax(ci, ifelse(is.na(sm), ci, as.integer(round(sm))))
  }
  dd2$.cls <- factor(labels[ci], levels = labels)
  g <- ggplot2::ggplot(dd2, ggplot2::aes(X_UTM33N, Y_UTM33N, fill = .cls)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::scale_fill_viridis_d(option = palette, direction = direction,
                                  drop = FALSE, na.translate = FALSE,
                                  guide = ggplot2::guide_legend(
                                    reverse = isTRUE(reverse_legend))) +
    ggplot2::labs(title = title, fill = legend, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  if (isTRUE(carto)) g <- g + do.call(carto_layers, c(list(dd), carto_args))
  g
}

# Classes of the q95 in cm (values range from 42 to 1000), with operational
# thresholds.
Q95_BREAKS <- c(0, 100, 200, 300, 500, Inf)
Q95_LABELS <- c("<100", "100-200", "200-300", "300-500", ">500")

# Light contours, without labels, for the uncertainty field of F43 (cm),
# inside the range of pi_width_smooth (about 39-471).
UNC_CONTOURS <- c(150, 250, 350)

# ---- Saving figures ---------------------------------------------------------
# The journal accepts vector graphics (PDF, EPS, with embedded fonts) or
# bitmaps of at least 300 dpi (2244 px at full page width) and 1000 dpi for
# line art. Therefore:
#  - charts (bars, curves, tree): vector PDF written with cairo_pdf, which
#    embeds the fonts and is the file LaTeX includes; next to it, the same
#    chart as a 300 dpi PNG for the Word copy, because pandoc does not embed
#    PDFs;
#  - maps (geom_raster): PNG at 600 dpi only. In a PDF the raster would enter
#    at the resolution of the grid, one cell per pixel, and some viewers blur
#    it.
DPI_CHART <- 300
DPI_MAP   <- 600
fig_file <- function(name, ext)
  file.path(PATHS$figs, paste0(sub("\\.(png|pdf)$", "", name), ".", ext))

# ggplot: vector PDF plus PNG; vector = FALSE (maps) writes the PNG only.
save_gg <- function(name, plot, w, h, vector = TRUE,
                    dpi = if (vector) DPI_CHART else DPI_MAP) {
  if (vector)
    ggplot2::ggsave(fig_file(name, "pdf"), plot, width = w, height = h,
                    device = grDevices::cairo_pdf)
  ggplot2::ggsave(fig_file(name, "png"), plot, width = w, height = h, dpi = dpi)
}

# Base graphics: draw() draws the figure and is called once per device, first
# for the vector PDF and then for the PNG.
save_base <- function(name, draw, w, h, dpi = DPI_CHART) {
  grDevices::cairo_pdf(fig_file(name, "pdf"), width = w, height = h)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(fig_file(name, "png"), width = w, height = h, units = "in",
                 res = dpi)
  tryCatch(draw(), finally = grDevices::dev.off())
}

save_map <- function(name, plot, w = 7, h = 6, dpi = DPI_MAP)
  save_gg(name, plot, w, h, vector = FALSE, dpi = dpi)

# Colours of the operating regimes (Figure 4), which the text calls green,
# amber and red. Green and amber are from Okabe and Ito; the red is dark. They
# stay distinguishable with deuteranopia and protanopia (minimum distance 11.4
# in OKLab x100, Machado et al. 2009 simulation) and, having different
# lightness, also when printed in grey.
COL_REGIME <- c(green = "#009E73", amber = "#E69F00", red = "#B2182B")

message("Setup OK | seed=", SEED, " | target=", TARGET)
