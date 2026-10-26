# 18_frontier_rules.R - Accuracy-fidelity frontier + explicit rules
# Produces (i) an explicit table of the frontier points and (ii) 3-4
# representative rules of the final e2tree for the main text.
source("R/00_setup.R")
suppressPackageStartupMessages({library(e2tree)})

# ---- 1. Accuracy vs nLoI frontier (explicit table) ---------------------------
# T15 is extended with the nLoI CIs over 10 subsamples when available.
fr <- read.csv(file.path(PATHS$tables, "T15_accuracy_fidelity_tradeoff.csv"),
               check.names = FALSE)
# For the selected configuration the CI is in T20:
ci <- read.csv(file.path(PATHS$tables, "T20_fidelity_ci.csv"),
               check.names = FALSE)
n_ci <- ci[ci$metric == "nLoI", ]
fr$nLoI_lo <- NA; fr$nLoI_hi <- NA
# 03c writes the selected configuration as "nt300 ns80 [selected]": without
# this row nLoI_lo and nLoI_hi of T40 would stay NA and the resampling
# interval would not reach the table of the paper.
sel <- grep("selected", fr$config)
if (length(sel) == 1 && nrow(n_ci) > 0) {
  fr$nLoI_lo[sel] <- round(n_ci$ci_low, 3)
  fr$nLoI_hi[sel] <- round(n_ci$ci_up, 3)
}
fr$RMSE_pct_loss <- round(100 * (fr$RF_RMSE - min(fr$RF_RMSE)) /
                          min(fr$RF_RMSE), 1)
fr <- fr[, c("config", "RF_RMSE", "RF_R2", "RMSE_pct_loss",
              "nLoI", "nLoI_lo", "nLoI_hi")]
write.csv(fr, file.path(PATHS$tables, "T40_frontier_explicit.csv"),
          row.names = FALSE)
cat("==== Accuracy / fidelity frontier ====\n"); print(fr, row.names = FALSE)

# ---- 2. Explicit rules of the final e2tree ----------------------------------
e2 <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))
tree <- e2$tree
tn <- e2tree::nodes(tree, terminal = TRUE)

# Extract the available columns
needed <- c("path", "n", "pred")
have <- intersect(needed, names(tn))
stopifnot(all(needed %in% have))

rules <- data.frame(
  rule_id = paste0("R", seq_len(nrow(tn))),
  n       = as.integer(as.character(tn$n)),
  pred_cm = round(as.numeric(as.character(tn$pred)), 1),
  path    = as.character(tn$path),
  stringsAsFactors = FALSE)
# Coarse semantic label (thin / medium / thick / very thick)
rules$class <- cut(rules$pred_cm, c(-Inf, 50, 100, 200, Inf),
                   labels = c("thin (<50)","medium (50-100)",
                              "thick (100-200)","very thick (>200)"))
rules <- rules[order(-rules$pred_cm), ]
write.csv(rules, file.path(PATHS$tables, "T41_e2tree_rules_clean.csv"),
          row.names = FALSE)
cat("\n==== E2Tree rules (ordered by predicted thickness) ====\n")
print(rules[, c("rule_id","n","pred_cm","class","path")], row.names = FALSE)

# Selection of the 4 representative examples: the thickest, two intermediate
# ones and the thinnest (largest n within each class).
sel_examples <- function(R) {
  pick <- integer(0)
  for (cl in levels(R$class)) {
    sub <- R[R$class == cl, , drop = FALSE]
    if (nrow(sub) == 0) next
    pick <- c(pick, which(R$rule_id == sub$rule_id[which.max(sub$n)]))
  }
  R[pick, ]
}
egs <- sel_examples(rules)
write.csv(egs, file.path(PATHS$tables, "T42_e2tree_rules_examples.csv"),
          row.names = FALSE)
cat("\n==== 4 representative rules ====\n")
print(egs[, c("rule_id","n","pred_cm","class","path")], row.names = FALSE)

saveRDS(list(frontier = fr, rules = rules, examples = egs),
        file.path(PATHS$models, "frontier_rules.rds"))
message("18_frontier_rules OK | rules=", nrow(rules))
