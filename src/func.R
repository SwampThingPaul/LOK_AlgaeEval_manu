
## FUNCTIONS
# st_txt now in AnalystHelper
## For NOAA CI data
ci.reverse.scaling.fun=function(DN){
  10^(3.0 / 250.0 * DN - 4.2)
}
ci.scaling.fun=function(ci){
  round(83.3 * (log10(ci[ci>0]) + 4.2))
}

model_fit_params2=function(Obs,Sim,DF=T){
  dat=na.omit(data.frame(Obs=Obs,Sim=Sim))
  Obs=dat$Obs
  Sim=dat$Sim
  MSE=1/N.obs(Obs)*sum((Obs-Sim)^2,na.rm=T)
  
  # res=Obs-Sim
  # fit.vals=Sim
  # rdf=N.obs(res)-2
  # mss=sum((fit.vals - mean(fit.vals,na.rm=T))^2,na.rm=T)
  # rss=sum(res^2,na.rm=T)
  # resvar=rss/rdf
  # sigma=sqrt(resvar);#RMSE
  # r2=mss/(mss+rss)
  # n <- sum(complete.cases(Sim))
  # r2 = 1 - (sum((Obs-Sim)^2, na.rm = T)/((n-1)*var(Obs, na.rm = T)));# "traditional" from caret package
  r2 = cor(Obs, Sim, use = "complete.obs")^2 ;# "corr" from caret package
  RMSE = sqrt(mean((Sim - Obs)^2, na.rm = T))
  MAE = mean(abs(Sim - Obs), na.rm = T)
  MAD = median(abs(Sim - Obs), na.rm = T)
  
  r.val=cov(Obs,Sim,use="pairwise.complete.obs")/(sd(Obs,na.rm=T)*sd(Sim,na.rm=T))
  # r2=r.val^2
  bias_mean=mean(Obs,na.rm=T)/mean(Sim,na.rm=T)
  sd_ratio=sd(Obs,na.rm=T)/sd(Sim,na.rm=T)
  MBE = mean(Sim - Obs,na.rm=T)
  # PBIAS=((mean(Sim,na.rm=T)-mean(Obs,na.rm=T))/mean(Sim,na.rm=T))*100
  
  KGE=1-sqrt((bias_mean-1)^2 + (sd_ratio-1)^2 + (r.val-1)^2)
  NSE = 1-MSE/sd(Obs)^2
  if(DF==T){rslt=data.frame(r.val=r.val,r2=r2,bias_mean=bias_mean,sd_ratio=sd_ratio,KGE=KGE,NSE=NSE,MAE = MAE,MAD=MAD,MBE=MBE,MSE=MSE,RMSE=RMSE)# ,PBIAS=PBIAS)
  }else{rslt=list(r.val=r.val,r2=r2,bias_mean=bias_mean,sd_ratio=sd_ratio,KGE=KGE,NSE=NSE,MAE = MAE,MSE=MSE,MAD=MAD,MBE=MBE,RMSE=RMSE)}#,PBIAS=PBIAS)}
  return(rslt)
}

dec.month=function(date){
  yr=as.numeric(format(date,"%Y"))
  leap=(yr%%4==0)&((yr%%100!=0)|(yr%%400 == 0))
  
  month_x <- month.abb[as.numeric(format(date,"%m"))]
  N_DAYS_IN_MONTHS <- c(
    Jan = 31L, Feb = 28L, Mar = 31L,
    Apr = 30L, May = 31L, Jun = 30L,
    Jul = 31L, Aug = 31L, Sep = 30L,
    Oct = 31L, Nov = 30L, Dec = 31L
  )
  
  n_days <- N_DAYS_IN_MONTHS[month_x]
  n_days[month_x == "Feb" & leap==T] <- 29L
  
  decmon=as.numeric(format(date,"%m")) + as.numeric(format(date,"%d"))/n_days
  return(decmon)
}

psum_fun <- function(...,na.rm=FALSE) { 
  rowSums(do.call(cbind,list(...)),na.rm=na.rm) } 

## for GAM periods of change
# download.file("https://gist.github.com/gavinsimpson/e73f011fdaaab4bb5a30/raw/82118ee30c9ef1254795d2ec6d356a664cc138ab/Deriv.R",
#               paste0(wd,"/src/Deriv.R"))
# source(paste0(wd,"/src/Deriv.R"))

## from .../Deriv.R
################################################
## Functions for derivatives of GAM(M) models ##
################################################
Deriv <- function(mod, n = 200, eps = 1e-7, newdata, term) {
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else { ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha/2), residual.df)
  ##for(i in term.labs[term]) {
  for(i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term,
                       eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  tVal <- qt(1 - (alpha/2), df.residual(x$gamModel))
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  CI <- confint(x, term = term)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  for(i in term) {
    upr <- CI[[i]]$upper
    lwr <- CI[[i]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}

dev.data.val=function(model,smooth.id,n=400,var.names,nc=nc){
  sm.term=smooths(model)
  
  df.res <- df.residual(model)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  
  newd <- gratia:::derivative_data(model,
                                   id = smooth.id, n = n,
                                   offset = NULL, order = 1L,
                                   type = "central", eps = 1e-7)
  terms.fit=predict(model,newd,type="terms",se.fit = T,nthreads=nc,discrete=T)
  tmp.fit=terms.fit$fit
  colnames(tmp.fit)=paste("fit",var.names,sep=".")
  tmp.SE=terms.fit$se.fit
  colnames(tmp.SE)=paste("SE",var.names,sep=".")
  newd=cbind(newd,tmp.fit,tmp.SE)
  newd$upper.CI=newd[,paste0("fit.",var.names[smooth.id])] + (crit.t* newd[,paste0("SE.",var.names[smooth.id])])
  newd$lower.CI=newd[,paste0("fit.",var.names[smooth.id])] - (crit.t* newd[,paste0("SE.",var.names[smooth.id])])
  
  model.d <- derivatives(model,newd,term=sm.term[smooth.id],type = "central",interval="confidence",ncores=nc)
  dsig <- signifD(newd[,paste0("fit.",var.names[smooth.id])],
                  d=model.d$derivative,# ver >0.8.1  d=model.d$.derivative,
                  upper=model.d$upper, # upper=model.d$.upper_ci,
                  lower=model.d$lower # lower=model.d$.lower_ci
  )
  dsig.incr=unlist(dsig$incr)
  dsig.decr=unlist(dsig$decr)
  newd=cbind(newd,dsig.incr,dsig.decr)
  return(newd)
}

stress.fit.fun=function(object){
  #source from vegan stressplot() function
  #https://github.com/vegandevs/vegan/blob/master/R/stressplot.R
  
  ## extract items to plot
  x <- object$diss
  y <- object$dist
  yf <- object$dhat
  
  ralscal <- 0
  ## Fit lines: linear (iregn=2) and hybrid (iregn=3) have a smooth line
  if (object$iregn > 1) {
    if (object$iregn == 3) {
      k <- seq(object$istart[2], object$ndis)
      yl <- range(yf[k])
      xl <- range(x[k])
      ralscal <- cor(y[k], yf[k])^2
    } else {
      yl <- range(yf)
      xl <- range(x)
      ralscal <- cor(y, yf)^2
    }
    #ln.val=data.frame(x=xl,y=yl);only for MonoMDS
  }
  ## Monotone line except in linear, and local has several...
  if (object$iregn != 2) {
    ist <- c(object$istart, object$ndis + 1)
    if (object$iregn == 3)
      object$ngrp <- 1
    for(j in 1:object$ngrp) {
      k <- seq(ist[j], ist[j+1]-1)
      ralscal <- ralscal + cor(y[k], yf[k])^2
      #ln.val=data.frame(x=x[k],y=yf[k])
    }
  }
  ## Stress as R2
  rstress <- 1 - object$stress^2
  ralscal <- if(object$iregn == 3) ralscal/2 else ralscal/object$ngrp
  Rst <- format(rstress, digits = 3)
  Ral <- format(ralscal, digits = 3)
  lab1 <- bquote("Non-metric fit, " * R^2 == .(Rst))
  lab2 <- bquote("Linear fit, " * R^2 == .(Ral))
  #return(list(data.frame(Fit=c("Non-metric","Linear"),R2=c(Rst,Ral)),ln.val))
  return(data.frame(Fit=c("Non-metric","Linear"),R2=c(Rst,Ral)))
}
glance2 <- function(x, ...) {
  s <- summary(x)
  
  data.frame(
    logLik = as.numeric(stats::logLik(x)),
    AIC = stats::AIC(x),
    BIC = stats::BIC(x),
    deviance = stats::deviance(x),
    df.residual = stats::df.residual(x),
    nobs = stats::nobs(x),
    npar = s$np,
    dispersion = s$scale,
    r.sq = s$r.sq,
    dev.exp = s$dev.expl
    
  )
}

LOOCV_fun <- function(model){
  h <- model$hat
  fitted_values <- predict(model, type = "response")
  y <- model$model[,1]
  
  # Compute residuals
  residuals <- residuals(model,type="response")
  
  # Approximate LOOCV predictions
  loocv_preds <- fitted_values - (residuals * h) / (1 - h)
  
  eval <- model_fit_params(y,loocv_preds)
  
  rslt <- list(
    pred = as.numeric(loocv_preds),
    actual = y,
    hat_diag = h,
    fit = eval
  )
  return(rslt)
}