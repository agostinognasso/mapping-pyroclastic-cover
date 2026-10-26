# 36_z0_sensitivity.R - 90% intervals of RF and the error of z0
# The QRF intervals are conditional on z0 as mapped: where the error of the
# isopach map is like that at the measured points, it is already in the spread
# from which QRF estimates the intervals. Far from the measurements the error
# can be larger, and there the coverage cannot be checked. This script runs a
# stress test: on the B15 folds of 17_pi_calibration.R, with the quantile
# configuration of the map (2000 trees, node size 20), the z0 of the held-out
# points receives an extra error, multiplicative and log-normal (sigma in
# log10 units), and the following are measured:
#  - plugin: the coverage of the current intervals, which do not know about
#    the error;
#  - propagated: the coverage of the mixture of the conditional distributions
#    over B draws of z0 around the observed value (Monte Carlo);
#  - propagated_noextra: the same without extra error, that is the cost of
#    propagating an error that the intervals already contain.
# Four error structures: independent point by point, a Gaussian field with a
# practical range of 5 and 20 km (exponential, on the distinct locations), and
# the same shift for the whole fold. The reference for sigma is Costa et al.
# (2012, par. 10): about 0.3 log10 units (a factor of 2) for the modelled
# isopach of the Campanian Ignimbrite, which here is a ceiling. Their fit
# excludes the measurements within 50 km of the eruptive centre.
# The R random stream used by quantreg stays that of 17 (the perturbations use
# a separate stream), so the base case reproduces Table 9: the script stops if
# it does not.
source("R/00_setup.R")
suppressPackageStartupMessages(library(ranger))

d  <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
pv <- predictor_vars(d)
f_cov <- as.formula(paste(TARGET, "~", paste(pv, collapse = "+")))
n_thr <- max(1L, parallel::detectCores() - 1L)
if (!file.exists(file.path(PATHS$models, "leave_massif.rds")))
  stop("Run 15_leave_massif.R first")
d$massif <- readRDS(file.path(PATHS$models, "leave_massif.rds"))$massif

set.seed(1)
b15 <- kmeans(scale(d[, c("X_UTM33N", "Y_UTM33N")]), 15, nstart = 5,
              iter.max = 50)$cluster

SIG    <- c(0.05, 0.10, 0.20, 0.30)       # standard deviation, log10 units
STRUCT <- c("iid", "grf5", "grf20", "fold")
R_REP  <- 10                              # replicates of the extra error
B_MC   <- 20                              # draws for the propagation
A      <- 0.10

winkler <- function(y, lo, hi) (hi - lo) +
  ifelse(y < lo, 2 / A * (lo - y), 0) + ifelse(y > hi, 2 / A * (y - hi), 0)
# Values of the terminal nodes (one draw per tree), as ranger does in
# predict(type = "quantiles"): rows = points, columns = trees.
nodevals <- function(m, nd) {
  tn <- predict(m, nd, type = "terminalNodes", num.threads = n_thr)$predictions + 1
  nv <- matrix(0, nrow(tn), ncol(tn))
  for (t in seq_len(ncol(tn))) nv[, t] <- m$random.node.values[tn[, t], t]
  nv
}
rowq <- function(M) t(apply(M, 1, quantile, c(0.05, 0.95), na.rm = TRUE))

# Error generator for the points of a fold. For the Gaussian fields the
# factorisation is done once per fold, on the distinct locations (coincident
# points have the same z0 and receive the same error).
make_draw <- function(te) {
  key <- paste(te$X_UTM33N, te$Y_UTM33N)
  u   <- !duplicated(key); id <- match(key, key[u])
  xy  <- as.matrix(te[u, c("X_UTM33N", "Y_UTM33N")])
  D   <- as.matrix(dist(xy))
  # With points a few metres apart and a 20 km range the matrix is nearly
  # singular: the minimum that suffices is added to the diagonal (at most 1% of
  # the variance, that is negligible independent noise).
  chol_j <- function(C) {
    for (j in c(1e-8, 1e-6, 1e-4, 1e-2)) {
      U <- tryCatch(chol(C + diag(j, nrow(C))), error = function(e) NULL)
      if (!is.null(U)) return(t(U))
    }
    stop("factorisation failed")
  }
  Lch <- lapply(c(grf5 = 5000, grf20 = 20000), function(L)
    chol_j(exp(-3 * D / L)))
  n <- nrow(te)
  function(st, s) switch(st,
    iid   = rnorm(n, 0, s),
    fold  = rep(rnorm(1, 0, s), n),
    s * as.vector(Lch[[st]] %*% rnorm(nrow(D)))[id])
}

acc <- list()
keep <- function(setting, s, st, r, te, lo, hi) {
  y <- te[[TARGET]]
  g <- data.frame(inside = y >= lo & y <= hi, width = hi - lo,
                  ws = winkler(y, lo, hi))
  ag <- rowsum(as.matrix(g), te$massif)
  acc[[length(acc) + 1]] <<- data.frame(setting = setting, sigma = s,
    struct = st, r = r, massif = as.integer(rownames(ag)),
    n = as.vector(table(te$massif)[rownames(ag)]), ag, row.names = NULL)
}

t0 <- Sys.time()
for (k in sort(unique(b15))) {
  trI <- which(b15 != k); teI <- which(b15 == k)
  m <- ranger(f_cov, d[trI, ], num.trees = 2000, mtry = 6, min.node.size = 20,
              quantreg = TRUE, seed = SEED, num.threads = n_thr)
  rs <- get(".Random.seed", envir = globalenv())
  set.seed(7000 + k)
  te <- d[teI, ]; n <- nrow(te)
  draw <- make_draw(te)
  q0 <- rowq(nodevals(m, te))
  keep("baseline", 0, "none", 0, te, q0[, 1], q0[, 2])
  for (s in SIG) for (st in STRUCT) for (r in 0:R_REP) {
    # r = 0: no extra error, only the propagation
    z0obs <- te$ThickIso * 10^(if (r == 0) 0 else draw(st, s))
    big <- te[rep(seq_len(n), B_MC + 1), ]
    big$ThickIso <- c(z0obs, unlist(lapply(seq_len(B_MC), function(b)
      z0obs * 10^draw(st, s))))
    nv <- nodevals(m, big)
    if (r > 0) {
      qp <- rowq(nv[seq_len(n), , drop = FALSE])
      keep("plugin", s, st, r, te, qp[, 1], qp[, 2])
    }
    pool <- do.call(cbind, lapply(seq_len(B_MC), function(b)
      nv[n * b + seq_len(n), , drop = FALSE]))
    qq <- rowq(pool)
    keep(if (r == 0) "propagated_noextra" else "propagated", s, st, r, te,
         qq[, 1], qq[, 2])
  }
  assign(".Random.seed", rs, envir = globalenv())
  cat(sprintf("fold %2d | n = %4d | %s\n", k, n,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1))))
}
X <- do.call(rbind, acc)

# ---- Summary ------------------------------------------------------------------
# Per replicate: coverage, width and Winkler aggregated over all folds; then
# mean (and range) over the replicates.
by_rep <- aggregate(cbind(n, inside, width, ws) ~ setting + sigma + struct + r,
                    X, sum)
by_rep$coverage <- by_rep$inside / by_rep$n
by_rep$width    <- by_rep$width / by_rep$n
by_rep$winkler  <- by_rep$ws / by_rep$n
S <- do.call(rbind, lapply(split(by_rep, with(by_rep,
       paste(setting, sigma, struct)), drop = TRUE), function(z)
  data.frame(setting = z$setting[1], sigma = z$sigma[1], struct = z$struct[1],
             n_rep = nrow(z), coverage = mean(z$coverage),
             coverage_min = min(z$coverage), coverage_max = max(z$coverage),
             width = mean(z$width), winkler = mean(z$winkler))))
S <- S[order(S$setting, S$struct, S$sigma), ]
base <- S[S$setting == "baseline", ]
S$width_ratio <- S$width / base$width

# The base case must reproduce Table 9 (T36)
t36 <- read.csv(file.path(PATHS$tables, "T36_pi_marginal.csv"))
cov36 <- t36$value[t36$metric == "coverage_90"]
wid36 <- t36$value[t36$metric == "mean_width_cm"]
if (abs(base$coverage - cov36) > 5e-4 || abs(base$width - wid36) > 0.05)
  stop(sprintf("Base case differs from T36: coverage %.4f against %.3f, width %.2f against %.2f",
               base$coverage, cov36, base$width, wid36))
cat(sprintf("\nBase case = T36: coverage %.3f, width %.1f cm, Winkler %.1f cm\n",
            base$coverage, base$width, base$winkler))

sect <- aggregate(cbind(n, inside) ~ setting + sigma + struct + massif, X, sum)
sect$coverage <- sect$inside / sect$n
sect$n <- sect$n / ifelse(sect$setting %in% c("baseline", "propagated_noextra"),
                          1, R_REP)
sect <- sect[order(sect$setting, sect$struct, sect$sigma, sect$massif),
             c("setting", "sigma", "struct", "massif", "n", "coverage")]

r3 <- function(z) { z[] <- lapply(z, function(v) if (is.numeric(v)) round(v, 4) else v); z }
write.csv(r3(S), file.path(PATHS$tables, "T53_z0_sensitivity.csv"), row.names = FALSE)
write.csv(r3(sect), file.path(PATHS$tables, "T53_z0_by_sector.csv"), row.names = FALSE)
saveRDS(list(summary = S, by_rep = by_rep, by_sector = sect,
             sigma = SIG, struct = STRUCT, R = R_REP, B = B_MC),
        file.path(PATHS$models, "z0_sensitivity.rds"))
cat("\n==== z0 stress test (90% coverage, mean over the replicates) ====\n")
print(r3(S[, c("setting", "sigma", "struct", "coverage", "coverage_min",
               "coverage_max", "width_ratio", "winkler")]), row.names = FALSE)

# ---- Figure: coverage against sigma -------------------------------------------
# Two series (current intervals, propagated error), one panel per structure.
# Okabe and Ito blue and vermilion, checked for deuteranopia and protanopia;
# the two series also have different shapes. No direct labels: with lines this
# close they would overlap, the legend is enough.
ST_LAB <- c(iid = "independent", grf5 = "correlated, 5 km",
            grf20 = "correlated, 20 km", fold = "same shift over the fold")
fd <- S[S$setting %in% c("plugin", "propagated"), ]
fd <- rbind(fd, transform(S[rep(which(S$setting == "baseline"), 2 * length(STRUCT)), ],
                          setting = rep(c("plugin", "propagated"), each = length(STRUCT)),
                          struct = rep(STRUCT, 2)))
fd$series <- factor(ifelse(fd$setting == "plugin", "current intervals",
                           "z0 error propagated"),
                    levels = c("current intervals", "z0 error propagated"))
fd$panel <- factor(ST_LAB[fd$struct], levels = ST_LAB)
COLS <- c("current intervals" = "#D55E00", "z0 error propagated" = "#0072B2")
SER_LAB <- c("current intervals", expression(z[0] * " error propagated"))
g <- ggplot(fd, aes(sigma, coverage, colour = series, shape = series)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.85, ymax = 0.95,
           fill = "grey92") +
  geom_hline(yintercept = 0.90, linetype = 2, colour = "grey35") +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_wrap(~ panel, nrow = 1) +
  scale_colour_manual(values = COLS, name = NULL, labels = SER_LAB) +
  scale_shape_manual(values = c(16, 17), name = NULL, labels = SER_LAB) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3)) +
  coord_cartesian(ylim = c(0.75, 1)) +
  labs(x = expression("Extra error of " * z[0] * ", standard deviation (log"[10] * " units)"),
       y = "Coverage of the 90% interval") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.spacing = grid::unit(1.4, "lines"))
save_gg("F54_z0_sensitivity", g, w = 9, h = 3.4)

# ---- LaTeX table for the supplement ------------------------------------------
f3 <- function(x) formatC(x, format = "f", digits = 3)
f1 <- function(x) formatC(x, format = "f", digits = 1)
pick <- function(st, s, set, col) S[[col]][S$struct == st & S$sigma == s & S$setting == set]
rows <- unlist(lapply(STRUCT, function(st) c(
  sprintf("\\multicolumn{8}{@{}l}{\\emph{%s}}\\\\", ST_LAB[[st]]),
  vapply(SIG, function(s) sprintf(
    "\\quad $%s$ & $\\times %s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$\\\\",
    formatC(s, format = "f", digits = 2), formatC(10^s, format = "f", digits = 2),
    f3(pick(st, s, "plugin", "coverage")), f3(pick(st, s, "propagated", "coverage")),
    f1(pick(st, s, "propagated", "width")), f1(pick(st, s, "plugin", "winkler")),
    f1(pick(st, s, "propagated", "winkler")),
    f3(pick(st, s, "propagated_noextra", "coverage"))), character(1)))))
writeLines(c(
  "% tab_supp_z0.tex -- generated by R/36_z0_sensitivity.R. Do not edit by hand.",
  sprintf("No extra error & -- & $%s$ & -- & $%s$ & $%s$ & -- & --\\\\",
          f3(base$coverage), f1(base$width), f1(base$winkler)),
  "\\midrule", rows),
  file.path(PATHS$paper, "tab_supp_z0.tex"))
message("36_z0_sensitivity OK")
