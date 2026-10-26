# 22_three_regimes.R - Three fit-for-purpose regimes
# The dichotomy "R2 > 0 = usable" is too coarse for operational use. Three
# bands are defined:
#   - direct operational use: R2 >= 0.40
#   - in-area support:        0.20 <= R2 < 0.40
#   - not usable:             R2 < 0.20
source("R/00_setup.R")
suppressPackageStartupMessages({library(ggplot2)})

SC <- readRDS(file.path(PATHS$models, "scenarios.rds"))
get_R2 <- function(sc, met = "rfsp") SC[[sc]]["R2", met]

scen <- data.frame(
  scenario = c("random5","sp40","sp15","sp5","leave_comune",
               "oper15_p10","oper15_p25","oper15_p50"),
  label = c("random 5-fold","spatial 40 blocks","spatial 15 blocks",
            "spatial 5 blocks","leave-municipality",
            "operational +10%","operational +25%","operational +50%"))
scen$R2_RFsp <- sapply(scen$scenario, get_R2)
scen$R2_GPT  <- sapply(scen$scenario, get_R2, met = "gpt")

# Classify by regime
classify <- function(r) {
  if (is.na(r)) return("unknown")
  if (r >= 0.40) return("operational use (>= 0.40)")
  if (r >= 0.20) return("in-area support (0.20-0.40)")
  return("not usable (< 0.20)")
}
scen$regime <- factor(sapply(scen$R2_RFsp, classify),
                       levels = c("operational use (>= 0.40)",
                                  "in-area support (0.20-0.40)",
                                  "not usable (< 0.20)"))

write.csv(scen, file.path(PATHS$tables, "T46_three_regimes.csv"),
          row.names = FALSE)
cat("==== Three operational regimes ====\n")
print(scen, row.names = FALSE)

# Figure F41 is drawn by R/13_figures_en.R, from the registry T47 of
# R/30_scenario_registry.R.

message("22_three_regimes OK")
