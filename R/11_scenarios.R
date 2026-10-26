# 11_scenarios.R - Matrix of validation scenarios x methods (honest CV)
# Produces output/models/scenarios.rds, used by the operational model (10).
source("R/00_setup.R")
suppressPackageStartupMessages({library(ranger);library(gstat);library(sp)})
d <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
# Covariates from predictor_vars() (00_setup.R), as in 01, 03 and 07: the 16
# variables of Table 2.
pv <- predictor_vars(d)
r2<-function(o,p) 1-sum((o-p)^2)/sum((o-mean(o))^2); rmse<-function(o,p) sqrt(mean((o-p)^2))
f_rf<-as.formula(paste("thick_cm~",paste(pv,collapse="+")))
f_sp<-as.formula(paste("thick_cm~",paste(c(pv,"X_UTM33N","Y_UTM33N"),collapse="+")))
# Same seed (SEED) as the rest of the pipeline, so a scenario has the same
# numbers in every table.
RF<-function(tr,f) ranger(f,tr,num.trees=400,mtry=6,min.node.size=5,seed=SEED,num.threads=4)
# Kriging (OK, RK, RFsp+RK). Precautions:
#  - residuals of points with coincident coordinates are aggregated: otherwise
#    the local kriging system is singular and krige() returns NA;
#  - no Gaussian model among the candidates: on these data it is numerically
#    unstable (predictions of the order of 1e15);
#  - the nugget starts from a fraction of the variance, which regularises the
#    system;
#  - empirical variogram with 250 m lags up to 10 km, as in
#    14_variogram_buffered_cv.R: with gstat's default lags the first one falls
#    beyond a kilometre and the short-range structure is lost;
#  - Sph and Exp from six starting values each; the fits that converge without
#    negative parameters are kept and the one with the smallest SSE is chosen;
#  - if no fit converges, or the predictions are not valid, the fold is left
#    without a kriging correction and the case is counted in KR_FAIL.
KR_WIDTH <- 250; KR_CUTOFF <- 10000
KR_FAIL <- new.env(); KR_FAIL$n <- 0L; KR_FAIL$why <- character(0)
KR_FAIL$fits <- list()   # one record per estimated variogram: what, model, range, nugget
fit_vgm <- function(v, s0) {
  best <- NULL
  for (m in c("Sph","Exp")) for (r0 in c(300, 1000, 3000)) for (nf in c(0.5, 0.9)) {
    w <- NULL
    f <- tryCatch(withCallingHandlers(
           fit.variogram(v, vgm(psill=(1-nf)*s0, model=m, range=r0, nugget=nf*s0),
                         warn.if.neg=FALSE),
           warning=function(e){ w <<- conditionMessage(e); invokeRestart("muffleWarning") }),
         error=function(e) NULL)
    ok <- !is.null(f) && is.null(w) && !isTRUE(attr(f,"singular")) &&
          all(f$psill>=0) && all(f$range>=0) && max(f$range)>0
    if (ok && (is.null(best) || attr(f,"SSErr") < attr(best,"SSErr"))) best <- f
  }
  best
}
kr<-function(tr,te,resid,what=""){
  bail <- function(why){ KR_FAIL$n <- KR_FAIL$n + 1L
                         KR_FAIL$why <- c(KR_FAIL$why, paste0(what, ": ", why))
                         rep(0,nrow(te)) }
  tryCatch({
    g <- aggregate(list(r=resid), by=list(x=tr$X_UTM33N, y=tr$Y_UTM33N), FUN=mean)
    coordinates(g)<-~x+y
    v  <- variogram(r~1, g, cutoff=KR_CUTOFF, width=KR_WIDTH)
    s0 <- stats::var(g$r, na.rm=TRUE)
    vm <- fit_vgm(v, s0)
    if (is.null(vm))                           return(bail("no converged fit"))
    st <- vm[vm$model!="Nug",]
    KR_FAIL$fits[[length(KR_FAIL$fits)+1]] <- data.frame(
      what=what, model=as.character(st$model),
      practical_range_m=round(if (st$model=="Exp") 3*st$range else st$range),
      nugget_share=round(vm$psill[vm$model=="Nug"]/sum(vm$psill),3))
    n<-data.frame(x=te$X_UTM33N,y=te$Y_UTM33N); coordinates(n)<-~x+y
    k<-suppressWarnings(krige(r~1,g,n,model=vm,nmax=40,debug.level=0)$var1.pred)
    if (all(is.na(k)))                         return(bail("all NA"))
    k <- ifelse(is.na(k),0,k)
    lim <- 3*max(abs(range(resid,na.rm=TRUE)))          # physical clamp
    if (max(abs(k))>lim){ k <- pmax(pmin(k,lim),-lim)
                          KR_FAIL$why <- c(KR_FAIL$why, paste0(what, ": clamped")) }
    k
  }, error=function(e) bail(paste("error:",conditionMessage(e)))) }
blocks<-function(K){set.seed(1);kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]),K,nstart=5,iter.max=50)$cluster}

eval_fold<-function(trI,teI){
  tr<-d[trI,]; te<-d[teI,]
  # Out-of-bag residuals for RK and RFsp+RK, as in the paper (Sections 3.1 and
  # 3.3): ranger returns them in $predictions. With in-bag predictions on the
  # training set the residuals would be shrunk towards zero.
  m<-RF(tr[,c("thick_cm",pv)],f_rf); rf_tr<-m$predictions; rf<-predict(m,te)$predictions
  ms<-RF(tr[,c("thick_cm",pv,"X_UTM33N","Y_UTM33N")],f_sp); rfsp<-predict(ms,te)$predictions
  rfsp_tr<-ms$predictions
  list(obs=te$thick_cm, gpt=d[teI,"GPT_model"],
       sapt=d[teI,"SAPT_model"], sept=d[teI,"SEPT_model"],
       ok = mean(tr$thick_cm)+kr(tr,te,tr$thick_cm-mean(tr$thick_cm),"ok"),
       rf = rf, rfsp = rfsp,
       rk = rf + kr(tr,te,tr$thick_cm-rf_tr,"rk"),
       rfsp_rk = rfsp + kr(tr,te,tr$thick_cm-rfsp_tr,"rfsp_rk"))
}
# The out-of-fold predictions are also kept point by point, in OOF (row of d,
# fold, training mean, one column per method), with the in-area points of the
# OP designs: 35_support_check.R uses them for the block means. The book-keeping
# does not use the random number generator, so T21 and T22 are unchanged.
OOF<-list()
run_cv<-function(folds,key){
  acc<-list(); idx<-integer(0); fk<-integer(0); nul<-numeric(0)
  for(k in sort(unique(folds))){ trI<-which(folds!=k); teI<-which(folds==k)
    r<-eval_fold(trI,teI)
    idx<-c(idx,teI); fk<-c(fk,rep(k,length(teI))); nul<-c(nul,rep(mean(d$thick_cm[trI]),length(teI)))
    for(nm in names(r)) acc[[nm]]<-c(acc[[nm]],r[[nm]]) }
  OOF[[key]]<<-list(pred=data.frame(i=idx,fold=fk,null=nul,acc),inarea=integer(0))
  meth<-setdiff(names(acc),"obs")
  sapply(meth,function(m) c(RMSE=rmse(acc$obs,acc[[m]]),R2=r2(acc$obs,acc[[m]])))
}
run_oper<-function(K,p,key){ # K blocks, but a fraction p of each block goes to training
  set.seed(1); b<-blocks(K); acc<-list(); idx<-integer(0); fk<-integer(0); nul<-numeric(0)
  ina<-integer(0)
  for(k in sort(unique(b))){ ink<-which(b==k); set.seed(k)
    add<-sample(ink,floor(p*length(ink))); teI<-setdiff(ink,add)
    if(length(teI)<5) next
    trI<-c(which(b!=k),add); r<-eval_fold(trI,teI)
    idx<-c(idx,teI); fk<-c(fk,rep(k,length(teI))); nul<-c(nul,rep(mean(d$thick_cm[trI]),length(teI)))
    ina<-c(ina,add)
    for(nm in names(r)) acc[[nm]]<-c(acc[[nm]],r[[nm]]) }
  OOF[[key]]<<-list(pred=data.frame(i=idx,fold=fk,null=nul,acc),inarea=ina)
  meth<-setdiff(names(acc),"obs")
  sapply(meth,function(m) c(RMSE=rmse(acc$obs,acc[[m]]),R2=r2(acc$obs,acc[[m]])))
}
set.seed(1); rnd<-sample(rep(1:5,length.out=nrow(d))); muni<-as.integer(factor(d$localitiy))
SC<-list(random5=run_cv(rnd,"random5"), sp40=run_cv(blocks(40),"sp40"),
  sp15=run_cv(blocks(15),"sp15"), sp5=run_cv(blocks(5),"sp5"),
  leave_comune=run_cv(muni,"leave_comune"),
  oper15_p10=run_oper(15,.10,"oper15_p10"), oper15_p25=run_oper(15,.25,"oper15_p25"),
  oper15_p50=run_oper(15,.50,"oper15_p50"))
saveRDS(OOF, file.path(PATHS$models, "oof_points.rds"))
cat("\n############ R2 by scenario x method ############\n")
R2<-t(sapply(SC,function(x) round(x["R2",],3))); print(R2)
cat("\n############ RMSE by scenario x method ############\n")
RM<-t(sapply(SC,function(x) round(x["RMSE",],1))); print(RM)
cat(sprintf("\nkriging: %d folds fell back (no correction)\n", KR_FAIL$n))
if (length(KR_FAIL$why)) print(table(KR_FAIL$why))
write.csv(as.data.frame(table(reason=KR_FAIL$why)),
          file.path(PATHS$tables, "T45_kriging_fallbacks.csv"), row.names=FALSE)
# The estimated variograms, one per fold and method, in a table from which
# the kriging range and nugget can be quoted.
KF <- do.call(rbind, KR_FAIL$fits)
write.csv(KF, file.path(PATHS$tables, "T45_kriging_fits.csv"), row.names=FALSE)
cat("\nkriging variograms (medians over folds):\n")
print(aggregate(cbind(practical_range_m, nugget_share) ~ what, KF, median))
print(table(KF$what, KF$model))
saveRDS(SC, file.path(PATHS$models, "scenarios.rds"))
write.csv(data.frame(scenario=rownames(R2),R2), file.path(PATHS$tables, "T21_scenarios_R2.csv"), row.names=FALSE)
write.csv(data.frame(scenario=rownames(RM),RM), file.path(PATHS$tables, "T22_scenarios_RMSE.csv"), row.names=FALSE)
