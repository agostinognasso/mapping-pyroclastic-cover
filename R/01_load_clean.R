# 01_load_clean.R - Read the data, clean it and apply the filters (as in
# Ebrahimi et al. 2024)
source("R/00_setup.R")

# ---- Reading ---------------------------------------------------------------
# sep = ";" with a decimal point: read.csv2 would use dec = "," and misread
# the numbers
raw <- read.csv(PATHS$data, sep = ";", dec = ".",
                check.names = FALSE, stringsAsFactors = FALSE)
n0 <- nrow(raw)

# Robust numeric coercion of covariates and response (handles 1.2E-06 notation)
num_cols <- setdiff(names(raw), c("note", "localitiy", "type"))
raw[num_cols] <- lapply(raw[num_cols], function(x) suppressWarnings(as.numeric(as.character(x))))

# Normalise the categorical fields
raw$note      <- tolower(trimws(raw$note))
raw$localitiy <- trimws(raw$localitiy)
raw$type      <- tolower(trimws(raw$type))

# ---- Filters, as in the paper ------------------------------------------------
# 1) Exclude the localities Forio and Procida (not on the mainland)
# 2) Exclude partial-thickness measurements (note == "parziale", the label
#    used in the published data)
excl_loc  <- raw$localitiy %in% c("Forio", "Procida")
excl_part <- raw$note == "parziale"

flt <- raw[!excl_loc & !excl_part, , drop = FALSE]

filter_log <- data.frame(
  step = c("dataset.csv, raw",
           "  excluded Forio/Procida",
           "  excluded partial measurements",
           "modelling dataset"),
  n = c(n0, n0 - sum(excl_loc), n0 - sum(excl_loc | excl_part), nrow(flt))
)

# ---- Missing values after coercion ------------------------------------------
na_after <- colSums(is.na(flt[num_cols]))
na_after <- na_after[na_after > 0]

# Rows with NA in the covariates or the response are dropped and logged
pv  <- predictor_vars(flt)
keep <- stats::complete.cases(flt[, c(TARGET, pv)])
n_dropped_na <- sum(!keep)
model_df <- flt[keep, , drop = FALSE]

filter_log <- rbind(filter_log,
  data.frame(step = "  rows with NA dropped (covariates/response)",
             n = nrow(model_df)))

# ---- Align with the prediction grid ------------------------------------------
newdata <- read.csv(PATHS$newdata, sep = ";", dec = ".",
                    check.names = FALSE, stringsAsFactors = FALSE)
newdata[] <- lapply(newdata, function(x) suppressWarnings(as.numeric(as.character(x))))
# Rename the grid columns to the training names:
#   Altitude  -> Elevation   (same variable, different name)
#   ThickIsop -> ThickIso    (z0, same variable, different name)
ren <- c(Altitude = "Elevation", ThickIsop = "ThickIso")
for (old in names(ren)) {
  new <- ren[[old]]
  if (old %in% names(newdata) && !new %in% names(newdata))
    names(newdata)[names(newdata) == old] <- new
}

pred_in_new <- setdiff(pv, names(newdata))   # covariates missing from the grid

# ---- Output -------------------------------------------------------------------
cat("\n==== FILTERS ====\n"); print(filter_log, row.names = FALSE)
cat("\nRF covariates (", length(pv), "):\n", sep = ""); print(pv)
cat("\nNA after numeric coercion (columns with >0):\n")
if (length(na_after)) print(na_after) else cat("  none\n")
cat("Rows dropped for NA:", n_dropped_na, "\n")
cat("\nCovariates missing from the prediction grid:",
    if (length(pred_in_new)) paste(pred_in_new, collapse = ", ") else "(none)", "\n")
cat("\nResponse thick_cm (modelling dataset):\n"); print(summary(model_df[[TARGET]]))

saveRDS(list(model_df = model_df, newdata = newdata, predictors = pv,
             filter_log = filter_log),
        file.path(PATHS$models, "data_clean.rds"))
write.csv(filter_log, file.path(PATHS$tables, "T0_filter_log.csv"), row.names = FALSE)
message("01_load_clean OK | n_model=", nrow(model_df), " | covariates=", length(pv))
