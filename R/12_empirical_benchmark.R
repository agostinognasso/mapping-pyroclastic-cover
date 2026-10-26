# 12_empirical_benchmark.R - Comparison of GPT, SAPT, SEPT and RFsp on the
# same spatial folds. The empirical models are columns of the dataset, with no
# training: the folds of 11_scenarios.R are rebuilt and the aggregate metrics
# are computed with the same logic (run_cv, run_oper).
source("R/00_setup.R")

d  <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
SC <- readRDS(file.path(PATHS$models, "scenarios.rds"))
r2 <- function(o,p) 1 - sum((o-p)^2)/sum((o-mean(o))^2)
rmse <- function(o,p) sqrt(mean((o-p)^2))
emp <- c(GPT="GPT_model", SAPT="SAPT_model", SEPT="SEPT_model")

blocks <- function(K){ set.seed(1)
  kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]), K, nstart=5, iter.max=50)$cluster }
pooled_part <- function(folds){      # partition CV: each point is tested once
  o<-d$thick_cm
  sapply(emp, function(col) c(RMSE=rmse(o,d[[col]]), R2=r2(o,d[[col]])))
}
pooled_oper <- function(K,p){        # as run_oper: test = (1-p) of each block
  set.seed(1); b<-blocks(K); idx<-c()
  for(k in sort(unique(b))){ ink<-which(b==k); set.seed(k)
    add<-sample(ink,floor(p*length(ink))); te<-setdiff(ink,add)
    if(length(te)>=5) idx<-c(idx,te) }
  o<-d$thick_cm[idx]
  sapply(emp, function(col) c(RMSE=rmse(o,d[[col]][idx]), R2=r2(o,d[[col]][idx])))
}

scen <- c("random5","sp40","sp15","sp5","leave_comune",
          "oper15_p10","oper15_p25","oper15_p50")
rows <- list()
for(s in scen){
  E <- if(grepl("^oper",s)) pooled_oper(15, as.numeric(sub(".*_p","",s))/100) else pooled_part()
  rows[[s]] <- data.frame(scenario=s,
    RFsp_R2  = round(SC[[s]]["R2","rfsp"],3),
    RFsp_RMSE= round(SC[[s]]["RMSE","rfsp"],1),
    GPT_R2 = round(E["R2","GPT"],3),  GPT_RMSE = round(E["RMSE","GPT"],1),
    SAPT_R2= round(E["R2","SAPT"],3), SAPT_RMSE= round(E["RMSE","SAPT"],1),
    SEPT_R2= round(E["R2","SEPT"],3), SEPT_RMSE= round(E["RMSE","SEPT"],1))
}
bench <- do.call(rbind, rows)
cat("==== Complete like-for-like benchmark (same folds) ====\n")
print(bench, row.names=FALSE)
write.csv(bench, file.path(PATHS$tables,"T28_fair_benchmark_full.csv"),
          row.names=FALSE)

# Figure: R2 by scenario, RFsp vs the three empirical models
library(ggplot2); library(tidyr)
L <- data.frame(scenario=factor(bench$scenario, levels=scen),
  RFsp=bench$RFsp_R2, GPT=bench$GPT_R2, SAPT=bench$SAPT_R2, SEPT=bench$SEPT_R2)
Lg <- pivot_longer(L, -scenario, names_to="model", values_to="R2")
p <- ggplot(Lg, aes(scenario, R2, fill=model)) +
  geom_col(position="dodge") + geom_hline(yintercept=0) +
  scale_fill_manual(values=c(RFsp="#1a9850",GPT="#d73027",
                             SAPT="#fc8d59",SEPT="#fee090")) +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=35,hjust=1)) +
  labs(title="Like-for-like benchmark: RFsp vs empirical models (same folds)",
       subtitle="The empirical models have R2 ~0 (or negative) in every regime",
       x=NULL, y=expression(R^2))
# Internal diagnostic. Figure 6 (F46_fair_benchmark_full) is written by
# 13_figures_en.R: a different name keeps a rerun of this script from
# overwriting it.
ggplot2::ggsave(file.path(PATHS$figs,"F46_fair_benchmark_full_diagnostic.png"), p,
                width=9, height=5, dpi=200)
message("12_empirical_benchmark OK")
