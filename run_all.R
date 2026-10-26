# run_all.R - Reproducible end-to-end pipeline
# Run from the repository root:  Rscript run_all.R
# Input data are in data/; results go to output/ and the LaTeX tables to
# output/tex/ (set PTCAM_OUT, PTCAM_PAPER and PTCAM_FIGS to write them
# elsewhere, see R/00_setup.R).
# The numbers 23 and 27 do not correspond to any script.
stopifnot(file.exists(file.path("data", "dataset.csv")))
steps <- c("R/01_load_clean.R", "R/02_eda.R", "R/03_rf_train.R",
           "R/03b_e2tree_setup.R", "R/04_e2tree.R", "R/04b_local_explain.R",
           "R/05_predict_maps.R", "R/09_quantile_uncertainty.R",
           "R/03c_tradeoff_nloi.R", "R/08_fidelity_ci.R",
           "R/07_spatial_cv.R", "R/11_scenarios.R",
           "R/12_empirical_benchmark.R", "R/10_operational_model.R",
           # Numbers for the Cilento quoted in the text, from the grid T27 and
           # from operational.rds produced by 10.
           "R/34_cilento_check.R",
           "R/06_report.R",
           # ---- Validation and comparison analyses --------------------
           "R/14_variogram_buffered_cv.R", "R/15_leave_massif.R",
           "R/16_empirical_recal.R", "R/17_pi_calibration.R",
           "R/18_frontier_rules.R", "R/19_z0_ablation.R",
           "R/20_ml_baselines.R", "R/21_measurement_sensitivity.R",
           "R/22_three_regimes.R",
           # Thresholds in cm and block means (reads oof_points.rds and T21
           # of 11), z0 stress test on the intervals (reads leave_massif.rds
           # of 15 and T36 of 17; about 35 minutes).
           "R/35_support_check.R", "R/36_z0_sensitivity.R",
           # Linear regression kriging on the same folds (reads T21 of 11,
           # which it must reproduce), sampled zone and area of
           # applicability (reads T27 of 10 and buffered_cv.rds of 14),
           # Figure 7 with the Wilson intervals (reads pi_calibration.rds of
           # 17). 30 and 13 read their outputs.
           "R/37_regression_kriging.R", "R/38_applicability.R",
           "R/39_coverage_figure.R",
           # From here the scripts read tables produced by the previous
           # steps; the figures (13) come after all the tables they use.
           # Surrogate CART and random partition: needs e2tree_fit.rds,
           # produced by 04.
           "R/32_cart_surrogate.R",
           # Absolute density of the local-data curve: before 13 (secondary
           # axis in pt/km^2) and 30 (columns of tab_curve). Reads
           # T27_operational_grid.csv, produced by 10.
           "R/31_curve_density.R",
           # Geometry of the CV designs (distance from the training set,
           # share within the range): reads T29_variogram_fit.csv, produced
           # by 14.
           "R/33_cv_geometry.R",
           # The registry runs before 13: Figures 4 and 6 (F41, F46) read
           # T47_scenario_registry.csv, to have the same acronyms and the
           # same order as the tables. It reads T44_curve_density (31) and
           # T47_cart_surrogate (32).
           "R/30_scenario_registry.R",
           "R/13_figures_en.R",
           # ---- Georeferencing and band reduction ---------------------
           "R/24_geotiff_export.R", "R/25_z0_debanding.R",
           "R/26_smooth_uncertainty.R",
           # ---- Catalogue of the e2tree rules (supplement) ------------
           "R/28_rule_catalogue.R",
           # ---- Supplement tables -------------------------------------
           "R/29_supp_tables.R")
for (s in steps) {
  cat("\n==========================================================\n")
  cat(">>> ", s, "  (", format(Sys.time()), ")\n", sep = "")
  cat("==========================================================\n")
  t <- system.time(source(s, local = new.env()))[3]
  cat(sprintf("<<< %s completed in %.1f s\n", s, t))
}
# Default folders: output/ and output/tex/ (R/00_setup.R)
cat("\nPipeline completed. Output in ",
    Sys.getenv("PTCAM_OUT", unset = "output"), "/, figures in ",
    Sys.getenv("PTCAM_FIGS", unset = "output/figs"), "/, LaTeX tables in ",
    Sys.getenv("PTCAM_PAPER", unset = "output/tex"), "/.\n", sep = "")
