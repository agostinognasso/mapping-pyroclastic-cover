# 09_quantile_uncertainty.R - Quantile regression forest: upper tail and uncertainty
# Forest only; e2tree plays no part here. It addresses the underestimation of
# extreme thicknesses, the most dangerous ones, and measures uncertainty with
# prediction intervals.
source("R/00_setup.R")

dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
fit <- readRDS(file.path(PATHS$models, "rf_fit.rds"))
su  <- readRDS(file.path(PATHS$models, "e2tree_setup.rds"))
cfg <- su$config
pv  <- fit$predictors
fml <- fit$formula
train <- fit$train[, c(TARGET, pv)]
test  <- fit$test[,  c(TARGET, pv)]

# Quantile RF (ranger): same family and settings as the explainable ensemble
qrf <- ranger(fml, data = train, num.trees = cfg$ntree,
               min.node.size = cfg$nodesize, mtry = cfg$mtry,
               quantreg = TRUE, seed = SEED,
               num.threads = max(1L, parallel::detectCores() - 1L))

qs  <- c(0.05, 0.5, 0.95)
pr  <- predict(qrf, test, type = "quantiles", quantiles = qs)$predictions
colnames(pr) <- c("q05", "q50", "q95")
y   <- test[[TARGET]]

# --- Coverage of the 90% prediction interval -----------------------------------
cover <- mean(y >= pr[, "q05"] & y <= pr[, "q95"])
piw   <- mean(pr[, "q95"] - pr[, "q05"])

# --- The dangerous tail: does q95 capture extreme thicknesses better? ----------
hi <- y >= quantile(y, 0.90)              # thickest deposits (risk)
rf_mean <- as.numeric(predict(fit$rf, test))
tail_tbl <- data.frame(
  metric = c("Tail MAE (RF mean)", "Tail MAE (RF q95)",
             "Tail bias, mean (pred-obs)", "Tail bias, q95 (pred-obs)",
             "PI90% coverage (all)", "Mean PI90% width (cm)"),
  value  = round(c(
    mae(y[hi], rf_mean[hi]), mae(y[hi], pr[hi, "q95"]),
    mean(rf_mean[hi] - y[hi]), mean(pr[hi, "q95"] - y[hi]),
    cover, piw), 2))
cat("==== DANGEROUS TAIL & UNCERTAINTY (ensemble only) ====\n")
print(tail_tbl, row.names = FALSE)
write.csv(tail_tbl, file.path(PATHS$tables, "T16_tail_uncertainty.csv"),
          row.names = FALSE)

# --- Maps: uncertainty and conservative (q95) map on the grid ------------------
grid <- dat$newdata
G <- grid[stats::complete.cases(grid[, pv]), ]
dom <- sapply(pv, function(v) range(fit$train[[v]], na.rm = TRUE))
G$in_domain <- apply(sapply(pv, function(v)
  G[[v]] >= dom[1, v] & G[[v]] <= dom[2, v]), 1, all)
gp <- predict(qrf, G, type = "quantiles", quantiles = qs)$predictions
G$q50 <- gp[, 1]; G$q95 <- gp[, 3]; G$pi_width <- gp[, 3] - gp[, 1]

theme_set(theme_minimal(base_size = 11))
Gin <- G[G$in_domain, ]
m1 <- ggplot(Gin, aes(X_UTM33N, Y_UTM33N, colour = pi_width)) +
  geom_point(size = 0.45) + coord_equal() +
  scale_colour_viridis_c(option = "inferno", direction = -1) +
  labs(title = "Ensemble uncertainty: width of the 90% interval (cm)",
       colour = "PI width", x = NULL, y = NULL) +
  theme(axis.text = element_blank())
m2 <- ggplot(Gin, aes(X_UTM33N, Y_UTM33N, colour = q95)) +
  geom_point(size = 0.45) + coord_equal() +
  scale_colour_viridis_c(option = "turbo") +
  labs(title = "Conservative scenario (0.95 quantile, cm) - in-domain cells only",
       colour = "q95 (cm)", x = NULL, y = NULL) +
  theme(axis.text = element_blank())
ggplot2::ggsave(file.path(PATHS$figs, "F26_map_uncertainty.png"), m1,
                width = 7, height = 6, dpi = 200)
ggplot2::ggsave(file.path(PATHS$figs, "F27_map_q95_conservative.png"), m2,
                width = 7, height = 6, dpi = 200)

saveRDS(list(qrf = qrf, tail = tail_tbl, qs = qs),
        file.path(PATHS$models, "qrf.rds"))
message("09_quantile_uncertainty OK | PI90% coverage=", round(cover, 3))
