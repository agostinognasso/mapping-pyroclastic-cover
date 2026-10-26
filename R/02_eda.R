# 02_eda.R - Exploratory analysis
source("R/00_setup.R")
dat <- readRDS(file.path(PATHS$models, "data_clean.rds"))
df  <- dat$model_df
pv  <- dat$predictors

theme_set(theme_minimal(base_size = 11))
# F01, F02 and F05, figures of the supplement, are written as vector PDF and
# as 300 dpi PNG (save_gg, 00_setup.R); F03 and F04, which are not in the
# paper, as 200 dpi PNG only.
ggsave2 <- function(name, p, w = 7, h = 5, vector = FALSE) {
  if (vector) save_gg(name, p, w, h)
  else ggplot2::ggsave(file.path(PATHS$figs, name), p, width = w, height = h,
                       dpi = 200)
}

# F01, F02 and F05 are supplement figures, which the journal publishes as
# received: no title in the figure (the caption carries it), covariates named
# as in Table 2 (the same names as in 29_supp_tables.R) and, in F05, blue and
# orange, which are distinguishable also by readers who cannot tell red from
# green.

# ---- 1. Distribution of the response -----------------------------------------
p1 <- ggplot(df, aes(thick_cm)) +
  geom_histogram(bins = 60, fill = "#2c7fb8", colour = "white") +
  labs(x = expression("Measured thickness" ~ italic(z) ~ "(cm)"), y = "Count")
p2 <- ggplot(df, aes(log1p(thick_cm))) +
  geom_histogram(bins = 60, fill = "#2c7fb8", colour = "white") +
  labs(x = expression(log(1 + italic(z))), y = "Count")
ggsave2("F01_target_hist.png", p1, vector = TRUE)
ggsave2("F02_target_log_hist.png", p2, vector = TRUE)

# ---- 2. Spatial distribution of the points -----------------------------------
p3 <- ggplot(df, aes(X_UTM33N, Y_UTM33N, colour = log1p(thick_cm))) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_viridis_c(option = "magma", direction = -1) +
  coord_equal() +
  labs(title = "Measurement points (colour = log1p thickness)",
       colour = "log1p(thick)")
ggsave2("F03_spatial_points.png", p3, w = 7, h = 6)

# ---- 3. Thickness by locality ------------------------------------------------
loc_ord <- df %>% group_by(localitiy) %>%
  summarise(med = median(thick_cm), n = n()) %>% arrange(desc(med))
df$localitiy <- factor(df$localitiy, levels = loc_ord$localitiy)
p4 <- ggplot(df, aes(localitiy, thick_cm)) +
  geom_boxplot(outlier.size = 0.4, fill = "grey85") +
  scale_y_log10() + coord_flip() +
  labs(title = "Thickness by municipality (log scale)", x = NULL, y = "thick_cm (cm)")
ggsave2("F04_thick_by_municipality.png", p4, w = 7, h = 6)

# ---- 4. Correlation of the covariates with the response ----------------------
cc <- sapply(pv, function(v) suppressWarnings(
  cor(df[[v]], df$thick_cm, use = "complete.obs", method = "spearman")))
cc <- sort(cc)
cor_df <- data.frame(var = factor(names(cc), levels = names(cc)), rho = cc)
# names of Table 2 of the paper, as plotmath (z0 and MSAVI2 with subscripts)
lab5 <- c(ThickIso = "italic(z)[0]", Elevation = "Altitude", Slope = "Slope",
          Aspect = "Aspect", Curvature = "Curvature",
          PlanCurv = "'Plan curvature'", ProfCurv = "'Profile curvature'",
          DistHydnet = "'Distance to hydrographic network'",
          FlowAcc = "'Flow accumulation'", FlowDir = "'Flow direction'",
          lnSPI = "ln(SPI)", STI = "STI", TWI = "TWI", NDVI = "NDVI",
          NCI = "NCI", MSAVI2 = "MSAVI[2]")
stopifnot(all(levels(cor_df$var) %in% names(lab5)))
p5 <- ggplot(cor_df, aes(rho, var, fill = rho > 0)) +
  geom_col() + scale_fill_manual(values = c("#d95f02", "#2c7fb8"), guide = "none") +
  scale_y_discrete(labels = function(v) parse(text = lab5[v])) +
  labs(x = expression("Spearman's" ~ rho ~ "with" ~ italic(z)), y = NULL)
ggsave2("F05_spearman_corr.png", p5, w = 7, h = 6, vector = TRUE)

# ---- 5. Descriptive tables ---------------------------------------------------
desc <- function(x) c(n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE),
                       sd = sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE),
                       median = median(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
desc_tab <- t(sapply(df[, c(TARGET, pv)], desc))
desc_tab <- round(as.data.frame(desc_tab), 3)
write.csv(desc_tab, file.path(PATHS$tables, "T1_descriptive.csv"))
write.csv(cor_df, file.path(PATHS$tables, "T2_spearman.csv"), row.names = FALSE)

cat("\n==== EDA ====\nDescriptives (response + covariates):\n")
print(desc_tab)
cat("\nLargest |rho| with the response:\n")
print(head(cor_df[order(-abs(cor_df$rho)), ], 8), row.names = FALSE)
message("02_eda OK | figures -> output/figs | tables -> output/tables")
