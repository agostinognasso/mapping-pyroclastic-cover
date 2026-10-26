# 04_e2tree.R - Explainable ensemble tree: post-hoc reading of the forest
#
# e2tree is not a predictor and does not replace the forest: it only serves to
# interpret it. Fidelity is measured with the structural measures of e2tree
# (nLoI, the normalised loss of information, and the eValidation measures),
# which compare the dissimilarity of the forest with the one rebuilt from the
# tree structure. Predictions of e2tree are never computed.
source("R/00_setup.R")
suppressPackageStartupMessages(library(future))

su  <- readRDS(file.path(PATHS$models, "e2tree_setup.rds"))
rf  <- su$rf_e2          # forest to interpret
S   <- su$S              # data used to induce the tree
pv  <- su$predictors
fml <- su$formula
setting <- su$setting

flatten_df <- function(df) {
  df <- as.data.frame(df)
  for (j in seq_along(df)) if (is.list(df[[j]]))
    df[[j]] <- vapply(df[[j]], function(z) paste(z, collapse = ";"), character(1))
  df
}

# ---- 1-2. Explanatory tree (cached: the induction is deterministic) ----------
cache <- file.path(PATHS$models, "e2tree_fit.rds")
prev  <- if (file.exists(cache)) readRDS(cache) else NULL
reuse <- !is.null(prev) && !is.null(prev$tree) &&
         identical(prev$setting, setting) &&
         identical(prev$config %||% NULL, su$config %||% NULL)
if (reuse) {
  tree <- prev$tree; D <- prev$D
  cat("e2tree: reusing the cached tree (setup unchanged)\n")
} else {
  t_dis <- system.time(
    D <- createDisMatrix(rf, data = S, label = TARGET,
                         parallel = list(active = FALSE, no_cores = 1)))[3]
  cat(sprintf("createDisMatrix: %.1f s  dim=%dx%d\n", t_dis, nrow(D), ncol(D)))
  t_tree <- system.time(tree <- e2tree(fml, S, D, rf, setting))[3]
  cat(sprintf("e2tree (level=%d): %.1f s\n", setting$level, t_tree))
}

# ---- 3. Structural fidelity (nLoI, Mantel and other measures) ---------------
# No predict(): the dissimilarity matrices of the forest and of the e2tree
# structure are compared. Low nLoI means high fidelity.
plan("sequential")                                  # no export to workers
options(future.globals.maxSize = 8 * 1024^3)
val <- tryCatch(
  eValidation(S, tree, D, test = "both", graph = FALSE,
              n_perm = 199, conf.level = 0.95, seed = SEED),
  error = function(e) { message("eValidation: ", conditionMessage(e)); NULL })

if (!is.null(val)) {
  meas <- val$measures
  cat("\n==== STRUCTURAL FIDELITY e2tree <-> ensemble (eValidation) ====\n")
  print(meas, row.names = FALSE)
  write.csv(flatten_df(meas), file.path(PATHS$tables, "T6_e2tree_fidelity.csv"),
            row.names = FALSE)
  nloi <- val$loi$nloi
  mant <- if (!is.null(val$mantel_test))
    c(z = val$mantel_test$z.stat, p = val$mantel_test$p) else NULL
  cat(sprintf("\nnLoI = %.4f  (0 = perfect fidelity; low = high fidelity)\n",
              nloi))
  if (!is.null(mant)) cat(sprintf("Mantel: z=%.3f  p=%.4f\n", mant[1], mant[2]))
  writeLines(capture.output(print(val)),
             file.path(PATHS$tables, "T8_e2tree_validation.txt"))
} else { meas <- NULL; nloi <- NA }

# ---- 4. Predictive accuracy, of the forest only ------------------------------
# e2tree is not a predictor, so it has no accuracy to measure.
fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
te  <- fit$test[, c(TARGET, pv)]
rf_te <- as.numeric(predict(rf, te))
ens_acc <- data.frame(model = "Random forest (ensemble, the only predictor)",
  set = "test", RMSE = round(rmse(te[[TARGET]], rf_te), 3),
  MAE = round(mae(te[[TARGET]], rf_te), 3),
  R2 = round(r2(te[[TARGET]], rf_te), 3))
cat("\n==== PREDICTIVE ACCURACY (RF ensemble only) ====\n")
print(ens_acc, row.names = FALSE)
write.csv(ens_acc, file.path(PATHS$tables, "T7_ensemble_accuracy.csv"),
          row.names = FALSE)

# ---- 5. Importance: internal reading of the ensemble through e2tree ----------
vi <- tryCatch(vimp(tree, data = S), error = function(e) NULL)
vi_e2 <- if (is.null(vi)) NULL else if (!is.null(vi$vimp)) vi$vimp else
         as.data.frame(vi[intersect(names(vi), c("Variable","MeanImpurityDecrease"))])
if (!is.null(vi_e2) && nrow(vi_e2)) {
  write.csv(flatten_df(vi_e2), file.path(PATHS$tables, "T9_e2tree_vimp.csv"),
            row.names = FALSE)
  cat("\n==== IMPORTANCE (e2tree, reading of the ensemble) ====\n")
  print(utils::head(vi_e2, 10))
} else message("vimp not available")

# ---- 6. Decision rules (terminal nodes = segments of the ensemble) -----------
tn  <- nodes(tree, terminal = TRUE)
write.csv(flatten_df(nodes(tree)), file.path(PATHS$tables, "T10_e2tree_nodes.csv"),
          row.names = FALSE)
write.csv(flatten_df(tn), file.path(PATHS$tables, "T11_e2tree_rules.csv"),
          row.names = FALSE)
cat("\nTerminal nodes (rules = segments of the ensemble's behaviour):",
    nrow(tn), "\n")
if (all(c("path", "pred") %in% names(tn))) {
  ord <- order(-as.numeric(tn$pred))
  cat("Rules (reference value of the segment, cm; a descriptive summary of",
      "the ensemble's behaviour, NOT a prediction):\n")
  for (i in ord) cat(sprintf("  ~%6.1f cm | n=%s | %s\n",
      as.numeric(tn$pred[i]), as.character(tn$n[i]), tn$path[i]))
}

# ---- 7. Figures ---------------------------------------------------------------
# Internal diagnostics only, with names different from those of the paper
# figures: F10_e2tree and F11_fidelity_measures are written by
# 13_figures_en.R, and the same name would overwrite them if this script were
# rerun after 13.
tryCatch({
  png(file.path(PATHS$figs, "F10_e2tree_diagnostic.png"), 1500, 1000, res = 150)
  plot_e2tree(tree, rf, main = "E2Tree - reading the random forest")
  dev.off()
}, error = function(e) message("plot_e2tree: ", conditionMessage(e)))
tryCatch(save_e2tree_html(tree, file.path(PATHS$figs, "F10_e2tree.html")),
         error = function(e) message("save_e2tree_html: ", conditionMessage(e)))

# Structural fidelity (eValidation measures), no comparison between predictions
if (!is.null(meas)) tryCatch({
  png(file.path(PATHS$figs, "F11_fidelity_measures_diagnostic.png"), 1100, 700, res = 150)
  op <- par(mar = c(8, 4.5, 4, 1))
  m <- meas[order(meas$observed), ]
  bp <- barplot(m$observed, names.arg = m$method, las = 2,
                col = "#2c7fb8", ylab = "observed value",
                main = "Structural fidelity e2tree <-> ensemble\n(native eValidation measures; low nLoI = high fidelity)")
  if (!is.null(m$ci_lower))
    arrows(bp, m$ci_lower, bp, m$ci_upper, angle = 90, code = 3,
           length = 0.04, col = "grey30")
  par(op); dev.off()
}, error = function(e) message("F11: ", conditionMessage(e)))

saveRDS(list(tree = tree, D = D, setting = setting, rf = rf,
             config = su$config %||% NULL,
             fidelity = meas, nloi = nloi, ensemble_accuracy = ens_acc,
             vimp = vi_e2, rules = tn),
        file.path(PATHS$models, "e2tree_fit.rds"))
message("04_e2tree OK | nLoI=", round(nloi, 4), " | nodes=", nrow(tn))
