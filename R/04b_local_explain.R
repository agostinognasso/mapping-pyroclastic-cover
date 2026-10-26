# 04b_local_explain.R - Local explanations: decision paths for selected cases
# The prediction comes from the forest. e2tree provides the explanation: which
# rule applies to the case and the reference value of that rule, which
# describes the forest's behaviour in that zone and is not an e2tree forecast.
source("R/00_setup.R")

fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
e2  <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))
tree <- e2$tree
rf_e2 <- e2$rf            # the interpreted forest, also used to predict
test <- fit$test
pv   <- fit$predictors

# Representative cases: high / median / low observed thickness ---------------
ord <- order(test[[TARGET]])
cases <- list(
  "low     (q05)"   = ord[max(1, floor(0.05 * length(ord)))],
  "median  (q50)"   = ord[floor(0.50 * length(ord))],
  "high    (q95)"   = ord[floor(0.95 * length(ord))],
  "extreme (max)"   = ord[length(ord)]
)

cat("==== LOCAL EXPLANATIONS (representative cases) ====\n")
asg  <- e2_assign(tree, test)                 # node/rule for each row
tnAll <- e2tree::nodes(tree, terminal = TRUE)
path_of <- function(nd) {
  j <- which(as.character(tnAll$node) == as.character(nd))
  if (length(j)) as.character(tnAll$path[j[1]]) else NA
}
out <- data.frame()
for (nm in names(cases)) {
  i  <- cases[[nm]]
  row <- test[i, , drop = FALSE]
  obs <- row[[TARGET]]
  rfp <- as.numeric(predict(rf_e2, row))      # prediction: the forest only
  node <- asg$node[i]                          # interpretation: which rule
  rule <- path_of(node)
  cat(sprintf("\n[%s]  observed=%.0f cm  |  RF ensemble PREDICTION=%.0f cm\n",
              nm, obs, rfp))
  cat(sprintf("  e2tree EXPLANATION: the ensemble's decision falls in segment/rule no. %s\n",
              as.character(node)))
  cat("  RULE:", if (is.na(rule)) "(root)" else rule, "\n")
  key <- c("ThickIso", "DistHydnet", "Elevation", "Slope", "NDVI", "NCI")
  key <- intersect(key, pv)
  print(round(row[, key, drop = FALSE], 2), row.names = FALSE)
  out <- rbind(out, data.frame(case = nm, observed = obs,
                 pred_ensemble_RF = round(rfp, 1),
                 segment_rule = as.character(node), rule = rule))
}
write.csv(out, file.path(PATHS$tables, "T13_local_explanations.csv"), row.names = FALSE)

# Explanation by municipality (if the localitiy column is in the test set) ----
if ("localitiy" %in% names(test)) {
  foc <- c("Cervinara", "Nocera Inferiore", "Vitulano")
  agg <- do.call(rbind, lapply(foc, function(L) {
    s <- test[test$localitiy == L, , drop = FALSE]
    if (!nrow(s)) return(NULL)
    a <- e2_assign(tree, s)
    dom_rule <- names(sort(table(a$node), decreasing = TRUE))[1]  # dominant rule
    data.frame(municipality = L, n = nrow(s),
               obs_med = median(s[[TARGET]]),
               pred_RF_med = round(median(as.numeric(predict(rf_e2, s))), 1),
               dominant_rule = dom_rule)
  }))
  cat("\n==== SUMMARY FOR SELECTED MUNICIPALITIES ====\n"); print(agg, row.names = FALSE)
  write.csv(agg, file.path(PATHS$tables, "T14_by_municipality.csv"), row.names = FALSE)
}
message("04b_local_explain OK")
