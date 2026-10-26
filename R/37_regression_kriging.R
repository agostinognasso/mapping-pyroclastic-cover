# 37_regression_kriging.R - Linear regression kriging on the same folds
# The comparison of Table 8 needs the classical geostatistical model with
# covariates: linear regression on the 16 covariates plus ordinary kriging of
# the residuals (regression kriging, Hengl et al. 2007; with the residuals of
# an OLS regression it is the operational form of kriging with external
# drift).
#
# Same designs, same folds and same kriging as 11_scenarios.R: the functions
# fit_vgm() and kr() are copied from there, unchanged. The check: OK
# recomputed here must reproduce the "ok" column of T21, otherwise the script
# stops. The regression residuals are those of the fit on the training set
# (OLS with 17 parameters on thousands of points: the shrinkage towards zero is
# negligible, unlike the in-bag residuals of a forest).
source("R/00_setup.R")
suppressPackageStartupMessages({library(gstat); library(sp)})
d  <- readRDS(file.path(PATHS$models, "data_clean.rds"))$model_df
pv <- predictor_vars(d)
f_lm <- as.formula(paste("thick_cm ~", paste(pv, collapse = " + ")))
r2s   <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
rmses <- function(o, p) sqrt(mean((o - p)^2))

# ---- Kriging: copy of 11_scenarios.R -------------------------------------------
KR_WIDTH <- 250; KR_CUTOFF <- 10000
KR_FAIL <- new.env(); KR_FAIL$n <- 0L; KR_FAIL$why <- character(0)
KR_FAIL$fits <- list()
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
    lim <- 3*max(abs(range(resid,na.rm=TRUE)))
    if (max(abs(k))>lim){ k <- pmax(pmin(k,lim),-lim)
                          KR_FAIL$why <- c(KR_FAIL$why, paste0(what, ": clamped")) }
    k
  }, error=function(e) bail(paste("error:",conditionMessage(e)))) }
blocks<-function(K){set.seed(1);kmeans(scale(d[,c("X_UTM33N","Y_UTM33N")]),K,nstart=5,iter.max=50)$cluster}

# ---- One fold: OLS, OLS + kriging of the residuals, OK as a check ------------
eval_fold <- function(trI, teI) {
  tr <- d[trI, ]; te <- d[teI, ]
  m  <- lm(f_lm, data = tr)
  lm_te <- predict(m, te)
  list(obs   = te$thick_cm,
       lm    = lm_te,
       lm_rk = lm_te + kr(tr, te, residuals(m), "lm_rk"),
       ok    = mean(tr$thick_cm) + kr(tr, te, tr$thick_cm - mean(tr$thick_cm), "ok"))
}
run_cv <- function(folds) {
  acc <- list()
  for (k in sort(unique(folds))) {
    r <- eval_fold(which(folds != k), which(folds == k))
    for (nm in names(r)) acc[[nm]] <- c(acc[[nm]], r[[nm]])
  }
  acc
}
run_oper <- function(K, p) {          # as in 11_scenarios.R
  set.seed(1); b <- blocks(K); acc <- list()
  for (k in sort(unique(b))) { ink <- which(b == k); set.seed(k)
    add <- sample(ink, floor(p * length(ink))); teI <- setdiff(ink, add)
    if (length(teI) < 5) next
    r <- eval_fold(c(which(b != k), add), teI)
    for (nm in names(r)) acc[[nm]] <- c(acc[[nm]], r[[nm]]) }
  acc
}
set.seed(1); rnd <- sample(rep(1:5, length.out = nrow(d)))
muni <- as.integer(factor(d$localitiy))
SC <- list(random5 = run_cv(rnd), sp40 = run_cv(blocks(40)),
           sp15 = run_cv(blocks(15)), sp5 = run_cv(blocks(5)),
           leave_comune = run_cv(muni),
           oper15_p10 = run_oper(15, .10), oper15_p25 = run_oper(15, .25),
           oper15_p50 = run_oper(15, .50))

OUT <- do.call(rbind, lapply(names(SC), function(s) {
  a <- SC[[s]]
  data.frame(scenario = s,
             R2_lm = r2s(a$obs, a$lm), R2_lm_rk = r2s(a$obs, a$lm_rk),
             R2_ok = r2s(a$obs, a$ok),
             RMSE_lm = rmses(a$obs, a$lm), RMSE_lm_rk = rmses(a$obs, a$lm_rk))
}))

# ---- Check: OK must reproduce T21 ---------------------------------------------
t21 <- read.csv(file.path(PATHS$tables, "T21_scenarios_R2.csv"))
dok <- abs(OUT$R2_ok - t21$ok[match(OUT$scenario, t21$scenario)])
if (any(dok > 0.0015))
  stop("OK does not reproduce T21 (folds or kriging differ from 11_scenarios.R): ",
       paste(OUT$scenario[dok > 0.0015], collapse = ", "))
cat("Check OK = T21 passed (maximum difference", signif(max(dok), 2), ")\n")

print(within(OUT, { R2_lm <- round(R2_lm, 3); R2_lm_rk <- round(R2_lm_rk, 3)
                    R2_ok <- round(R2_ok, 3); RMSE_lm <- round(RMSE_lm, 1)
                    RMSE_lm_rk <- round(RMSE_lm_rk, 1) }), row.names = FALSE)
KF <- do.call(rbind, KR_FAIL$fits)
cat(sprintf("\nkriging: %d folds without correction\n", KR_FAIL$n))
if (length(KR_FAIL$why)) print(table(KR_FAIL$why))
cat("\nvariograms of the OLS residuals (medians over folds):\n")
print(aggregate(cbind(practical_range_m, nugget_share) ~ what, KF, median))

write.csv(OUT, file.path(PATHS$tables, "T55_regression_kriging.csv"), row.names = FALSE)
write.csv(KF, file.path(PATHS$tables, "T55_regression_kriging_fits.csv"), row.names = FALSE)
write.csv(as.data.frame(table(reason = KR_FAIL$why)),
          file.path(PATHS$tables, "T55_regression_kriging_fallbacks.csv"), row.names = FALSE)
saveRDS(list(summary = OUT, fits = KF, fail = KR_FAIL$n),
        file.path(PATHS$models, "regression_kriging.rds"))
message("37_regression_kriging OK")
