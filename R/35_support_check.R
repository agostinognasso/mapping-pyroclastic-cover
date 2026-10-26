# 35_support_check.R - Thresholds in centimetres and errors on block means
#
# A. The regime thresholds (R^2 0.40 and 0.20, Section 4.1) translated into
#    RMSE: RMSE = SD * sqrt(1 - R^2), with the SD of the points evaluated in
#    each design. Next to it, the errors in cm of RF and RFsp for each design.
# B. Section 5.3 says that users need volumes per catchment, whereas the
#    validation is on points. Here the out-of-fold predictions of
#    11_scenarios.R (oof_points.rds) are averaged over square blocks of 0.5,
#    1, 2 and 5 km, and the observed and predicted block means are compared.
#    Averaging cancels the independent part of the error, not the systematic
#    part: the share that remains is MSE of the means / MSE of the points. The
#    measures are in cm and not in R^2, because averaging reduces the variance
#    and the R^2 of the means is not comparable with that of the points.
#    In the designs with in-area data (OP) the mean of the in-area points of
#    the same block is also compared, the simplest estimator a user holding a
#    survey would use.
# Limits, also stated in the supplement: the "truth" of a block is the mean at
# the measured points, not the mean over the area; in the OP designs the
# in-area points are a random subsample of the same clusters, which favours
# the local mean; the blocks all lie in surveyed zones.
source("R/00_setup.R")

d   <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
OOF <- readRDS(file.path(PATHS$models, "oof_points.rds"))
z   <- d[[TARGET]]

# ---- Check: the saved predictions reproduce T21 ------------------------------
t21 <- read.csv(file.path(PATHS$tables, "T21_scenarios_R2.csv"))
for (k in names(OOF)) for (m in c("rf", "rfsp", "ok")) {
  p <- OOF[[k]]$pred
  v <- round(r2(p$obs, p[[m]]), 3); w <- t21[[m]][t21$scenario == k]
  if (abs(v - w) > 0.0015) stop("oof_points.rds does not reproduce T21: ", k, " ", m)
}
stopifnot(all(sapply(OOF, function(o) all(o$pred$obs == z[o$pred$i]))))

# ---- A. Thresholds in centimetres -------------------------------------------
SD_ALL <- sd(z)
thr <- data.frame(R2 = c(0.40, 0.20, 0),
                  label = c("operational use", "in-area support", "no skill"))
thr$RMSE_cm <- SD_ALL * sqrt(1 - thr$R2)
thr$RMSE_reduction_pct <- 100 * (1 - sqrt(1 - thr$R2))
cat(sprintf("SD of the thickness (n = %d): %.2f cm\n", length(z), SD_ALL))
print(thr, row.names = FALSE)

DES <- c(random5 = "R5", sp5 = "B5", sp15 = "B15", sp40 = "B40",
         leave_comune = "LOMO", oper15_p10 = "OP10", oper15_p25 = "OP25",
         oper15_p50 = "OP50")
cm_tab <- do.call(rbind, lapply(names(DES), function(k) {
  p <- OOF[[k]]$pred
  sdk <- sqrt(mean((p$obs - mean(p$obs))^2))      # SD (population) of the evaluated points
  do.call(rbind, lapply(c(RF = "rf", RFsp = "rfsp"), function(m) {
    e <- p[[m]] - p$obs
    data.frame(design = DES[[k]], key = k, model = names(which(c(RF = "rf", RFsp = "rfsp") == m)),
               n = nrow(p), R2 = r2(p$obs, p[[m]]), RMSE = sqrt(mean(e^2)),
               MAE = mean(abs(e)), MedAE = median(abs(e)),
               within50_pct = 100 * mean(abs(e) <= 50),
               SD_scored = sdk,
               RMSE_at_040 = sdk * sqrt(0.60), RMSE_at_020 = sdk * sqrt(0.80))
  }))
}))
# The constant prediction (training mean of the fold) as a reference
cm_null <- do.call(rbind, lapply(names(DES), function(k) {
  p <- OOF[[k]]$pred; e <- p$null - p$obs
  data.frame(design = DES[[k]], MAE_null = mean(abs(e)), MedAE_null = median(abs(e)))
}))
cm_tab <- merge(cm_tab, cm_null, by = "design", sort = FALSE)
write.csv(cbind(cm_tab[, 1:4], round(cm_tab[, -(1:4)], 3)),
          file.path(PATHS$tables, "T51_threshold_cm.csv"), row.names = FALSE)
write.csv(transform(thr, RMSE_cm = round(RMSE_cm, 2),
                    RMSE_reduction_pct = round(RMSE_reduction_pct, 2),
                    SD_cm = round(SD_ALL, 2)),
          file.path(PATHS$tables, "T51_threshold_cm_equivalents.csv"), row.names = FALSE)
cat("\n==== Errors in cm by design ====\n")
print(cbind(cm_tab[, 1:3], round(cm_tab[, c("R2", "RMSE", "MAE", "MedAE",
       "within50_pct", "SD_scored", "RMSE_at_040", "MAE_null")], 2)), row.names = FALSE)

# ---- B. Block means -----------------------------------------------------------
L_SET   <- c(500, 1000, 2000, 5000)          # block side (m)
MIN_PTS <- 10                                # evaluated points per block
MIN_LOC <- 3                                 # in-area points for the local mean
N_BOOT  <- 1000
# Order of the tables of the paper; B5 and LOMO show what happens to blocks
# far from the training data.
BDES    <- c(sp5 = "B5", sp15 = "B15", sp40 = "B40", leave_comune = "LOMO",
             oper15_p10 = "OP10", oper15_p25 = "OP25", oper15_p50 = "OP50")
MODS    <- c(RF = "rf", RFsp = "rfsp", OK = "ok", "RF+RK" = "rk")
uid_of  <- function(L) paste(floor(d$X_UTM33N / L), floor(d$Y_UTM33N / L))

# For a design and a side length: one row per block with the means and, for
# each model, the sum of squared errors at the points (for the uncancelled
# share) and the local mean of the in-area points.
block_table <- function(k, L, min_pts = MIN_PTS) {
  p <- OOF[[k]]$pred; uid <- uid_of(L); p$u <- uid[p$i]
  p <- p[p$u %in% names(which(table(p$u) >= min_pts)), ]
  U <- data.frame(u = sort(unique(p$u)))
  U$n   <- as.vector(table(p$u)[U$u])
  U$obs <- as.vector(tapply(p$obs, p$u, mean)[U$u])
  U$null <- as.vector(tapply(p$null, p$u, mean)[U$u])
  for (m in MODS) {
    U[[m]] <- as.vector(tapply(p[[m]], p$u, mean)[U$u])
    U[[paste0("sse_", m)]] <- as.vector(tapply((p[[m]] - p$obs)^2, p$u, sum)[U$u])
  }
  ia <- OOF[[k]]$inarea
  if (length(ia)) {
    U$loc  <- as.vector(tapply(z[ia], uid[ia], mean)[U$u])
    U$nloc <- as.vector(table(uid[ia])[U$u]); U$nloc[is.na(U$nloc)] <- 0
  } else { U$loc <- NA_real_; U$nloc <- 0 }
  U
}

stats_of <- function(U, m) {
  e <- U[[m]] - U$obs
  mse_u <- mean(e^2); mse_p <- sum(U[[paste0("sse_", m)]]) / sum(U$n)
  c(RMSE_pt = sqrt(mse_p), RMSE_unit = sqrt(mse_u),
    rel_err_pct = 100 * median(abs(e) / U$obs),
    unexplained_pct = 100 * mse_u / mse_p,
    skill_vs_null = 1 - mse_u / mean((U$null - U$obs)^2),
    R2_unit = r2(U$obs, U[[m]]),
    slope = unname(coef(lm(U$obs ~ U[[m]]))[2]))
}
# Simple bootstrap of the blocks. Stratifying it by municipality does not
# work: with 5 km blocks many municipalities have a single block, which
# resampling never varies, and the intervals come out narrow for no reason.
boot_idx <- function(U) sample.int(nrow(U), nrow(U), replace = TRUE)
ci <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)

set.seed(SEED)
res <- list(); loc <- list(); pair <- list()
for (k in names(BDES)) for (L in L_SET) {
  U <- block_table(k, L)
  B <- replicate(N_BOOT, boot_idx(U), simplify = FALSE)
  for (mn in names(MODS)) {
    m <- MODS[[mn]]; st <- stats_of(U, m)
    bs <- t(vapply(B, function(ix) stats_of(U[ix, ], m)[c("RMSE_unit", "rel_err_pct",
                                                          "unexplained_pct")], numeric(3)))
    res[[length(res) + 1]] <- data.frame(design = BDES[[k]], key = k, L_m = L,
      model = mn, n_units = nrow(U), n_pts = sum(U$n), t(st),
      RMSE_unit_lo = ci(bs[, 1])[1], RMSE_unit_hi = ci(bs[, 1])[2],
      rel_err_lo = ci(bs[, 2])[1], rel_err_hi = ci(bs[, 2])[2],
      unexpl_lo = ci(bs[, 3])[1], unexpl_hi = ci(bs[, 3])[2])
  }
  # RF against RFsp on the same units (difference of the RMSE of the means)
  rm_u <- function(ix, m) sqrt(mean((U[[m]][ix] - U$obs[ix])^2))
  dd <- vapply(B, function(ix) rm_u(ix, "rf") - rm_u(ix, "rfsp"), numeric(1))
  pair[[length(pair) + 1]] <- data.frame(design = BDES[[k]], L_m = L, n_units = nrow(U),
    RMSE_RF = rm_u(seq_len(nrow(U)), "rf"), RMSE_RFsp = rm_u(seq_len(nrow(U)), "rfsp"),
    diff_RF_minus_RFsp = rm_u(seq_len(nrow(U)), "rf") - rm_u(seq_len(nrow(U)), "rfsp"),
    diff_lo = ci(dd)[1], diff_hi = ci(dd)[2])
  # Local mean, on the blocks with at least MIN_LOC in-area points
  V <- U[!is.na(U$loc) & U$nloc >= MIN_LOC, ]
  if (nrow(V) >= 5) {
    BV <- replicate(N_BOOT, boot_idx(V), simplify = FALSE)
    rv <- function(ix, m) sqrt(mean((V[[m]][ix] - V$obs[ix])^2))
    row <- data.frame(design = BDES[[k]], L_m = L, n_units = nrow(V),
                      median_inarea_pts = median(V$nloc),
                      RMSE_local = rv(seq_len(nrow(V)), "loc"))
    for (mn in names(MODS)) row[[paste0("RMSE_", mn)]] <- rv(seq_len(nrow(V)), MODS[[mn]])
    dl <- vapply(BV, function(ix) rv(ix, "rf") - rv(ix, "loc"), numeric(1))
    row$diff_RF_minus_local <- row$RMSE_RF - row$RMSE_local
    row$diff_lo <- ci(dl)[1]; row$diff_hi <- ci(dl)[2]
    row$rel_err_local_pct <- 100 * median(abs(V$loc - V$obs) / V$obs)
    row$rel_err_RF_pct    <- 100 * median(abs(V$rf - V$obs) / V$obs)
    loc[[length(loc) + 1]] <- row
  }
}
RES  <- do.call(rbind, res);  rownames(RES) <- NULL
PAIR <- do.call(rbind, pair); LOC <- do.call(rbind, loc)

# Sensitivity to the minimum number of points per block (1 km, RF)
sens <- do.call(rbind, lapply(names(BDES), function(k) do.call(rbind,
  lapply(c(5, 10, 20), function(mp) { U <- block_table(k, 1000, mp)
    data.frame(design = BDES[[k]], min_pts = mp, n_units = nrow(U),
               t(stats_of(U, "rf")[c("RMSE_unit", "rel_err_pct", "unexplained_pct")])) }))))

rnd <- function(x) { x[] <- lapply(x, function(v) if (is.numeric(v)) round(v, 3) else v); x }
write.csv(rnd(RES),  file.path(PATHS$tables, "T52_support_check.csv"), row.names = FALSE)
write.csv(rnd(PAIR), file.path(PATHS$tables, "T52_support_rf_vs_rfsp.csv"), row.names = FALSE)
write.csv(rnd(LOC),  file.path(PATHS$tables, "T52_support_localmean.csv"), row.names = FALSE)
write.csv(rnd(sens), file.path(PATHS$tables, "T52_support_minpts.csv"), row.names = FALSE)
saveRDS(list(thresholds = thr, cm = cm_tab, support = RES, rf_vs_rfsp = PAIR,
             localmean = LOC, minpts = sens, sd = SD_ALL),
        file.path(PATHS$models, "support_check.rds"))

cat("\n==== Block means, 1 km ====\n")
print(rnd(RES[RES$L_m == 1000, c("design", "model", "n_units", "RMSE_pt", "RMSE_unit",
  "RMSE_unit_lo", "RMSE_unit_hi", "rel_err_pct", "unexplained_pct", "slope")]), row.names = FALSE)
cat("\n==== RF against RFsp on block means ====\n"); print(rnd(PAIR), row.names = FALSE)
cat("\n==== Mean of the in-area points against the models ====\n"); print(rnd(LOC), row.names = FALSE)
cat("\n==== Sensitivity to the minimum points per block (RF, 1 km) ====\n"); print(rnd(sens), row.names = FALSE)

# ---- Figure: error of the means against the support --------------------------
# x axis: from the point to 5 km blocks; one curve per model, one panel per
# design. Okabe and Ito colours (blue, vermilion, green, amber), checked for
# the most common forms of colour blindness; amber has little contrast on
# white, so the four series also have different shapes, and the values are in
# the tables of the supplement.
SUP <- c("point", "0.5 km", "1 km", "2 km", "5 km")
pt <- do.call(rbind, lapply(names(BDES), function(k) do.call(rbind,
  lapply(names(MODS), function(mn) { p <- OOF[[k]]$pred
    data.frame(design = BDES[[k]], model = mn, support = "point",
               RMSE = sqrt(mean((p[[MODS[[mn]]]] - p$obs)^2)), lo = NA, hi = NA) }))))
bl <- with(RES, data.frame(design, model, support = paste(L_m / 1000, "km"),
                           RMSE = RMSE_unit, lo = RMSE_unit_lo, hi = RMSE_unit_hi))
FIG_DES <- c("B5", "B15", "OP10", "OP25", "OP50")   # B40 and LOMO in the table only
fd <- rbind(pt, bl)
fd <- fd[fd$design %in% FIG_DES, ]
fd$support <- factor(fd$support, levels = SUP)
fd$design  <- factor(fd$design, levels = FIG_DES)
fd$model   <- factor(fd$model, levels = names(MODS))
COLS <- c(RF = "#0072B2", RFsp = "#D55E00", OK = "#009E73", "RF+RK" = "#E69F00")
g <- ggplot(fd, aes(support, RMSE, colour = model, shape = model, group = model)) +
  geom_line(linewidth = 0.6, position = position_dodge(width = 0.35)) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5,
                 position = position_dodge(width = 0.35), na.rm = TRUE) +
  geom_point(size = 2, position = position_dodge(width = 0.35)) +
  facet_wrap(~ design, nrow = 1) +
  scale_colour_manual(values = COLS, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Support of the prediction (measurement point or square block)",
       y = "RMSE (cm)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_gg("F53_support_curve", g, w = 9, h = 3.6)

# ---- LaTeX tables for the supplement ------------------------------------------
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
# Design and number of blocks only on the first row of the group. No
# \multicolumn at the top of the file: after \input it would no longer be the
# start of the cell and LaTeX would stop ("Misplaced \omit").
r1 <- RES[RES$L_m == 1000, ]
lines1 <- unlist(lapply(seq_along(BDES), function(j) {
  ds <- BDES[[j]]; s <- r1[r1$design == ds, ]
  c(if (j > 1) "\\addlinespace",
    sprintf("%s & %s & %s & $%s$ & $%s$ \\,[$%s$, $%s$] & $%s$ & $%s$ & $%s$ & $%s$\\\\",
            c(ds, rep("", nrow(s) - 1)), c(sprintf("$%d$", s$n_units[1]), rep("", nrow(s) - 1)),
            s$model, f1(s$RMSE_pt), f1(s$RMSE_unit), f1(s$RMSE_unit_lo), f1(s$RMSE_unit_hi),
            f1(s$rel_err_pct), f1(s$unexplained_pct), f2(s$skill_vs_null), f2(s$slope)))
}))
writeLines(c("% tab_supp_support.tex -- generated by R/35_support_check.R. Do not edit by hand.",
             lines1), file.path(PATHS$paper, "tab_supp_support.tex"))
# Errors in cm by design (part A), RF and RFsp side by side
cmw <- merge(cm_tab[cm_tab$model == "RF", ], cm_tab[cm_tab$model == "RFsp", ],
             by = "design", suffixes = c("_RF", "_RFsp"), sort = FALSE)
cmw <- cmw[match(c("R5", "B5", "B15", "B40", "LOMO", "OP10", "OP25", "OP50"), cmw$design), ]
writeLines(c("% tab_supp_cm.tex -- generated by R/35_support_check.R. Do not edit by hand.",
  sprintf("%s & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$\\\\", cmw$design,
          f1(cmw$RMSE_RF), f1(cmw$MAE_RF), f1(cmw$MedAE_RF),
          f1(cmw$RMSE_RFsp), f1(cmw$MAE_RFsp), f1(cmw$MedAE_RFsp), f1(cmw$MAE_null_RF))),
  file.path(PATHS$paper, "tab_supp_cm.tex"))
l1 <- LOC[LOC$L_m == 1000, ]
writeLines(c("% tab_supp_localmean.tex -- generated by R/35_support_check.R. Do not edit by hand.",
  sprintf("%s & $%d$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ \\,[$%s$, $%s$]\\\\",
          l1$design, l1$n_units, formatC(round(l1$median_inarea_pts), format = "d"),
          f1(l1$RMSE_local), f1(l1$RMSE_RF), f1(l1[["RMSE_RF+RK"]]), f1(l1$RMSE_RFsp),
          f1(l1$RMSE_OK), f1(l1$diff_RF_minus_local), f1(l1$diff_lo), f1(l1$diff_hi))),
  file.path(PATHS$paper, "tab_supp_localmean.tex"))
message("35_support_check OK")
