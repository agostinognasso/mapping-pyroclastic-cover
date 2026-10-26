# 06_report.R - Assemble the final summary (tables + figure index)
source("R/00_setup.R")

rd <- function(f) { p <- file.path(PATHS$tables, f)
  if (file.exists(p)) read.csv(p, check.names = FALSE) else NULL }

cat("================  SUMMARY OF RESULTS  ================\n\n")

cat("## Dataset filters\n");          print(rd("T0_filter_log.csv"), row.names = FALSE)
cat("\n## Test-set performance (RF vs empirical models vs published)\n")
print(rd("T4_test_performance.csv"), row.names = FALSE)
cat("\n## INTERPRETIVE fidelity: e2tree reconstructs the ensemble\n")
print(rd("T6_e2tree_fidelity.csv"), row.names = FALSE)
cat("\n## PREDICTIVE accuracy: RF ensemble ONLY (e2tree does not predict)\n")
print(rd("T7_ensemble_accuracy.csv"), row.names = FALSE)
cat("\n## RF ensemble importance (top)\n")
vir <- rd("T5_rf_importance.csv"); if (!is.null(vir)) print(utils::head(vir, 8), row.names = FALSE)
cat("\n## e2tree importance - internal reading of the ensemble (top)\n")
vie <- rd("T9_e2tree_vimp.csv");  if (!is.null(vie)) print(utils::head(vie, 8), row.names = FALSE)

# ---- Accuracy <-> fidelity trade-off figure ------------------------------------
to <- rd("T15_accuracy_fidelity_tradeoff.csv")
if (!is.null(to) && "nLoI" %in% names(to)) {
  d <- to
  d$isP  <- grepl("predictive", d$config)
  d$isCh <- grepl("selected", d$config)
  d$col  <- ifelse(d$isP, "#d73027", ifelse(d$isCh, "#1a9850", "#2c7fb8"))
  d$pchv <- ifelse(d$isP, 4, ifelse(d$isCh, 17, 19))
  ymax <- max(d$nLoI, na.rm = TRUE) * 1.15

  # Internal diagnostic: F30_tradeoff, the supplement figure, is written by
  # 13_figures_en.R.
  png(file.path(PATHS$figs, "F30_tradeoff_diagnostic.png"), 1200, 800, res = 150)
  op <- par(mar = c(5, 5.5, 4, 2), xpd = NA)
  # Inverted y axis: low nLoI (a more faithful tree) at the top
  plot(d$RF_RMSE, d$nLoI, type = "n",
       xlim = c(min(d$RF_RMSE) - 0.6, max(d$RF_RMSE) + 0.8),
       ylim = c(ymax, 0),
       xlab = "Forest accuracy: test RMSE (cm), worse towards the right",
       ylab = "nLoI (top: more faithful and interpretable)",
       main = "Frontier between forest accuracy and e2tree fidelity (nLoI)")
  fr <- d[order(d$RF_RMSE), ]
  lines(fr$RF_RMSE, fr$nLoI, col = "grey60", lty = 2)
  points(d$RF_RMSE, d$nLoI, pch = d$pchv, col = d$col,
         cex = ifelse(d$isCh, 2.4, 1.8), lwd = 2)
  text(d$RF_RMSE, d$nLoI, d$config, pos = 3, cex = 0.72, col = d$col,
       offset = 0.9)
  legend("bottomright", bty = "n", cex = 0.8,
         pch = c(4, 19, 17), col = c("#d73027", "#2c7fb8", "#1a9850"),
         legend = c("Predictive RF (high nLoI = low fidelity)",
                    "Regularised RFs", "Selected RF"))
  par(op); dev.off()
  cat("\n## Accuracy (ensemble) <-> nLoI fidelity trade-off\n")
  print(to, row.names = FALSE)
}

figs <- sort(list.files(PATHS$figs, pattern = "\\.png$"))
cat("\n## Figures written (", length(figs), ")\n", sep = "")
for (f in figs) cat("  -", f, "\n")

# sessionInfo, for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(PATHS$paper, "sessionInfo.txt"))
cat("\nsessionInfo -> sessionInfo.txt in the LaTeX output folder\n")
message("06_report OK")
