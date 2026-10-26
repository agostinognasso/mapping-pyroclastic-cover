# 28_rule_catalogue.R - Catalogue of the rules of the explanatory tree (e2tree)
#
# The zone map (F45) colours each cell with the node of the tree, but on its
# own it does not say what "rule 64" is: which thresholds define it and which
# thickness it indicates. The script extracts the terminal nodes, simplifies
# the path from the root to the leaf into readable intervals (merges the
# conditions on the same variable and writes "!X<=t" as "X>t") and produces:
#   - PATHS$tables/T37_rule_catalogue.csv  (plain text, for the pipeline)
#   - PATHS$paper/tab_rules.tex            (LaTeX rows for the supplement)
# ThickIso is the initial isopach thickness z0 (notation of the paper).
source("R/00_setup.R")

e2  <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))
tn  <- e2tree::nodes(e2$tree, terminal = TRUE)
T27 <- read.csv(file.path(PATHS$tables, "T27_operational_grid.csv"),
                check.names = FALSE)
cov <- table(factor(T27$rule_zone))
covpct <- 100 * as.numeric(cov) / sum(cov); names(covpct) <- names(cov)

rnd <- function(x) if (abs(x) >= 10) as.character(round(x)) else as.character(round(x, 2))

# Raw path -> vector of simplified conditions (one interval per variable).
simplify_conds <- function(p) {
  parts <- trimws(strsplit(p, " & ", fixed = TRUE)[[1]])
  ub <- list(); lb <- list(); ordv <- c()
  for (cond in parts) {
    neg <- startsWith(cond, "!"); c2 <- sub("^!", "", cond)
    m <- regmatches(c2, regexec("^([[:alnum:]_]+)[[:space:]]*<=[[:space:]]*(-?[0-9.]+)$", c2))[[1]]
    if (length(m) == 0) next
    v <- m[2]; val <- as.numeric(m[3]); ordv <- c(ordv, v)
    if (!neg) ub[[v]] <- min(c(ub[[v]], val)) else lb[[v]] <- max(c(lb[[v]], val))
  }
  vars <- unique(ordv)
  sapply(vars, function(v) {
    hi <- ub[[v]]; lo <- lb[[v]]
    if (!is.null(hi) && !is.null(lo)) sprintf("%s < %s <= %s", rnd(lo), v, rnd(hi))
    else if (!is.null(hi)) sprintf("%s <= %s", v, rnd(hi))
    else sprintf("%s > %s", v, rnd(lo))
  }, USE.NAMES = FALSE)
}

# Text condition in LaTeX: variables in \mathrm{}, "<=" becomes "\le".
# Column names stay those of the dataset; the paper prints the symbols of its
# notation: ThickIso is z0, the initial thickness, not z.
DISPLAY_NAMES <- c(ThickIso = "z_{0}", Elevation = "Altitude")
cond_to_latex <- function(cond) {
  s <- gsub("([A-Za-z][A-Za-z0-9_]*)", "\\\\mathrm{\\1}", cond)
  s <- gsub("<=", "\\\\le", s)
  # \mathrm{ThickIso} -> z_{0} (a maths symbol, not \mathrm)
  s <- gsub("\\\\mathrm\\{ThickIso\\}", "z_{0}", s)
  s <- gsub("\\\\mathrm\\{Elevation\\}", "\\\\mathrm{Altitude}", s)
  paste0("$", s, "$")
}

conds <- lapply(as.character(tn$path), simplify_conds)
df <- data.frame(
  node     = tn$node,
  cm       = round(tn$pred),
  n        = tn$n,
  grid_pct = as.numeric(covpct[as.character(tn$node)]),
  rule     = vapply(conds, paste, "", collapse = ", "),
  rule_tex = vapply(conds, function(cc) paste(vapply(cc, cond_to_latex, ""),
                                              collapse = ", "), ""),
  stringsAsFactors = FALSE
)
df$grid_pct[is.na(df$grid_pct)] <- 0
df <- df[order(-df$cm), ]

# ---- CSV (pipeline) ---------------------------------------------------------
write.csv(data.frame(node = df$node, expected_thickness_cm = df$cm,
                     n_train = df$n, grid_pct = round(df$grid_pct, 1),
                     rule = df$rule),
          file.path(PATHS$tables, "T37_rule_catalogue.csv"), row.names = FALSE)

# ---- LaTeX fragment (\input in the paper) ---------------------------------
gp <- ifelse(df$grid_pct < 0.1, "$<0.1$",
             formatC(df$grid_pct, format = "f", digits = 1))
rows <- sprintf("%d & %d & %d & %s & %s \\\\",
                df$node, df$cm, df$n, gp, df$rule_tex)
writeLines(rows, file.path(PATHS$paper, "tab_rules.tex"))

cat("Rule catalogue: ", nrow(df), " nodes | thicknesses ",
    min(df$cm), "-", max(df$cm), " cm\n", sep = "")
message("28_rule_catalogue OK | T37_rule_catalogue.csv + paper/tab_rules.tex")
