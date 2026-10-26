# 08_fidelity_ci.R - Robustness of the fidelity: nLoI over K subsamples
# Structural fidelity (nLoI) depends on the subsample used to induce the tree.
# This script characterises its variability: mean and 95% CI over K replicates
# (fixed explainable configuration; only the subsample seed changes).
source("R/00_setup.R")
suppressPackageStartupMessages(library(future)); plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

su  <- readRDS(file.path(PATHS$models, "e2tree_setup.rds"))
fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
cfg <- su$config; pv <- su$predictors; fml <- su$formula
train <- fit$train[, c(TARGET, pv)]

KREP  <- as.integer(Sys.getenv("E2_KREP", "10"))
SUBN  <- as.integer(Sys.getenv("E2_CI_SUBN", "1500"))

strat_idx <- function(y, n, seed) {
  set.seed(seed)
  br <- unique(quantile(y, seq(0, 1, 0.1), na.rm = TRUE))
  q  <- cut(y, br, include.lowest = TRUE)
  unlist(lapply(split(seq_along(y), q), function(ix)
    sample(ix, max(1, round(n * length(ix) / length(y))))), use.names = FALSE)
}

vals <- data.frame()
for (r in seq_len(KREP)) {
  sd_seed <- SEED + r
  S <- train[strat_idx(train[[TARGET]], SUBN, sd_seed), ]
  set.seed(sd_seed)
  rfo <- randomForest(fml, S, ntree = cfg$ntree, nodesize = cfg$nodesize,
                      mtry = cfg$mtry)
  Dm <- createDisMatrix(rfo, data = S, label = TARGET,
                        parallel = list(active = FALSE, no_cores = 1))
  t2 <- e2tree(fml, S, Dm, rfo,
               list(impTotal = 0.1, maxDec = 1e-6, n = 2, level = cfg$level))
  ev <- eValidation(S, t2, Dm, test = "both", graph = FALSE,
                     n_perm = 99, seed = sd_seed)
  mant <- if (!is.null(ev$mantel_test)) ev$mantel_test$z.stat else NA
  vals <- rbind(vals, data.frame(rep = r, nLoI = ev$loi$nloi,
                                 mantel_z = mant,
                                 term_nodes = nrow(nodes(t2, terminal = TRUE))))
  cat(sprintf("rep %2d/%d | nLoI=%.4f | mantel z=%.3f | nodes=%d\n",
              r, KREP, vals$nLoI[r], vals$mantel_z[r], vals$term_nodes[r]))
}

ci <- function(x) { x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x); n <- length(x)
  c(mean = m, sd = s, ci_low = m - 1.96 * s / sqrt(n),
    ci_up = m + 1.96 * s / sqrt(n), min = min(x), max = max(x)) }
summ <- rbind(nLoI = ci(vals$nLoI), mantel_z = ci(vals$mantel_z),
              term_nodes = ci(vals$term_nodes))
summ <- round(summ, 4)
cat("\n==== nLoI FIDELITY - robustness over", KREP, "subsamples ====\n")
print(summ)
write.csv(vals, file.path(PATHS$tables, "T19_fidelity_reps.csv"), row.names = FALSE)
write.csv(data.frame(metric = rownames(summ), summ),
          file.path(PATHS$tables, "T20_fidelity_ci.csv"), row.names = FALSE)
message("08_fidelity_ci OK | mean nLoI=", round(mean(vals$nLoI, na.rm=TRUE), 4))
