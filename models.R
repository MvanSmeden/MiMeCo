################################################################
######################### MODELS ###############################
################################################################

###############################################
# Author: M van Smeden                        #
# Date first version: March 27 2019           #
# Simulations confounding, measurement error  #
# and missing values                          #
###############################################

################## catch and summary functions ############################ 
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){
    W <<- w
    invokeRestart("muffleWarning")}
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler),warning = W)
}

lm.summ <- function(f,v,lab,conflevel){
  CI <- confint(f$value,level=conflevel)[v,]
  est <- f$value$coefficients[v]
    TEMP <- data.frame(warn =ifelse(is.null(f$warning),"no",f$warning),
               est=est,
               err=est-arg$beta,
               sqr.err= (est-arg$beta)^2,
               estse=vcov(f$value)[v,v],
               ci.covering=dplyr::between(arg$beta,as.numeric(CI[1]),as.numeric(CI[2])),
               ci.width= as.numeric(CI[2])-as.numeric(CI[1]),
               ci.lb= as.numeric(CI[1]),
               ci.ub= as.numeric(CI[2]),
               N = nrow(model.frame(f$value)),
               R2 = summary(f$value)$r.squared
    )
    rownames(TEMP) <- NULL
    colnames(TEMP) <- paste0(lab,".",colnames(TEMP))
  TEMP
}

mice.summ <- function(fs,v,lab,err,conflevel){
  pooled <-(pool(fs))
  est <- pooled$pooled[v,]
  CI <- as.numeric(summary(pooled, conf.int = TRUE,conf.level=conflevel)[v,6:7])
    TEMP <- data.frame(warn =ifelse(is.null(err),"no",unlist(err)),
                     est=est[1,"estimate"],
                     err=est[1,"estimate"]-arg$beta,
                     sqr.err=(est[1,"estimate"]-arg$beta)^2,
                     estse=summary(pooled)[v,"std.error"],
                     ci.covering=dplyr::between(arg$beta,CI[1],CI[2]),
                     ci.width= CI[2]-CI[1],
                     ci.lb= as.numeric(CI[1]),
                     ci.ub= as.numeric(CI[2]),
                     N=pooled$pooled[1,"dfcom"]+arg$P+2,
                     R2=pool.r.squared(fs)[1,"est"],
                     nimp = length(fs$analyses)
    )
    rownames(TEMP) <- NULL
    colnames(TEMP) <- paste0(lab,".",colnames(TEMP))
  TEMP
}

lavaan.summ <- function(f,v,lab,conflevel){
  coeff <- parameterEstimates(f$value,level=conflevel)
    coeff.id <- which(coeff$lhs=="Y"&coeff$rhs==v)
     est <- as.numeric(coeff[coeff.id,"est"])
     CI <- as.numeric(coeff[coeff.id,c("ci.lower","ci.upper")])
      TEMP <- data.frame(warn =ifelse(is.null(f$warning),"no",unlist(f$warning$message)),
                        est=est,
                        err=est-arg$beta,
                        sqr.err=(est-arg$beta)^2,
                        estse= as.numeric(coeff[coeff.id,"se"]),
                        ci.covering=dplyr::between(arg$beta,CI[1],CI[2]),
                        ci.width= CI[2]-CI[1],
                        ci.lb= as.numeric(CI[1]),
                        ci.ub= as.numeric(CI[2]),
                        N=inspect(f$value,'ntotal'),
                        R2=inspect(f$value,'r2')["Y"]
       )
       rownames(TEMP) <- NULL
       colnames(TEMP) <- paste0(lab,".",colnames(TEMP))
    TEMP
}

mecor.summ <- function(f,v,lab,conflevel){
  est <- as.numeric(f$value$corfit[[1]][v])
  se <- sqrt(as.numeric(f$value$corvar[[1]][v]))
  z <- qnorm(1-(1-conflevel)/2)
  CI <- c(as.numeric(est-z*se),as.numeric(est+z*se))
  TEMP <- data.frame(warn =ifelse(is.null(f$warning),"no",unlist(f$warning$message)),
                     est=est,
                     err=est-arg$beta,
                     sqr.err=(est-arg$beta)^2,
                     estse= se,
                     ci.covering=dplyr::between(arg$beta,CI[1],CI[2]),
                     ci.width= CI[2]-CI[1],
                     ci.lb= as.numeric(CI[1]),
                     ci.ub= as.numeric(CI[2]))
  rownames(TEMP) <- NULL
  colnames(TEMP) <- paste0(lab,".",colnames(TEMP))
  TEMP
}

imp.summ <- function(OUT,f,imp,lab,conflevel){
  est <- mean(OUT[,"est"])
  se <-  sqrt(mean(OUT[,"estse"]^2)+var(OUT[,"est"])+var(OUT[,"est"])/nrow(OUT))
  z <- qnorm(1-(1-conflevel)/2)
  CI <- c(as.numeric(est-z*se),as.numeric(est+z*se))
  TEMP <- data.frame(warnMI =ifelse(is.null(imp$warning),"no",unlist(imp$warning$message)),
                     warnRC=ifelse(any(sapply(f,function(x)!is.null(x$warning))),sum(sapply(f,function(x)!is.null(x$warning))),"no"),
                     est=est,
                     err=est-arg$beta,
                     sqr.err=(est-arg$beta)^2,
                     estse= se,
                     ci.covering=dplyr::between(arg$beta,CI[1],CI[2]),
                     ci.width= CI[2]-CI[1],
                     ci.lb= as.numeric(CI[1]),
                     ci.ub= as.numeric(CI[2]))
  rownames(TEMP) <- NULL
  colnames(TEMP) <- paste0(lab,".",colnames(TEMP))
  TEMP
}

jags.summ <- function(f,v,lab,conflevel,nchains){
  mcsum <- summary(f$value,quantiles=c((1-conflevel)/2,1-(1-conflevel)/2))
  MCsamplesV <- as.numeric(sapply(f$value,function(x)x[,v]))
  dens <- density(MCsamplesV)
  est <- dens$x[which.max(dens$y)]
  CI <- quantile(MCsamplesV,c((1-conflevel)/2,1-(1-conflevel)/2))
  TEMP <- data.frame(warn =ifelse(is.null(f$warning),"no",unlist(f$warning$message)),
                     est=est,
                     estMean= mean(MCsamplesV),
                     estMedian= median(MCsamplesV),
                     err=est-arg$beta,
                     sqr.err=(est-arg$beta)^2,
                     estse= sd(MCsamplesV), #mcsum$statistics[v,'Naive SE'], # updated 5/7/2019 by Bas
                     ci.covering=dplyr::between(arg$beta,CI[1],CI[2]),
                     ci.width= CI[2]-CI[1],
                     ci.lb= as.numeric(CI[1]),
                     ci.ub= as.numeric(CI[2]))
  if(nchains>1){
    GRstat <-gelman.diag(f$value)
    TEMP$GRstat_A <- GRstat$psrf[v,1]
    TEMP$GRstat_mv <- GRstat$mpsrf
  } 
  rownames(TEMP) <- NULL
  colnames(TEMP) <- paste0(lab,".",colnames(TEMP))
  TEMP
}

################## models ################################## 
gs.ols <- function(D,arg){
  set.seed(arg$seed.dgm)
  Dgs <- simdgm(arg)
  f <- tryCatch.W.E(lm(Y~.,data=Dgs))
  data.frame(arg,lm.summ(f,v="A",lab="gs",conflevel=arg$conflevel))
}

ccunivar.ols <- function(D,arg){
  Dcomp <- D[complete.cases(D),]
  f <- tryCatch.W.E(lm(Y~Astar,data=Dcomp))
  data.frame(arg,lm.summ(f,v="Astar",lab="ccunivar",conflevel=arg$conflevel))
}

univar.ols <- function(D,arg){
  f <- tryCatch.W.E(lm(Y~Astar,data=D))
  data.frame(arg,lm.summ(f,v="Astar",lab="univar",conflevel=arg$conflevel))
}

naive.ols <- function(D,arg){
   f <- tryCatch.W.E(lm(paste0('Y~Astar+',paste(paste0('L.',1:arg$P),collapse='+')),data=D))
    data.frame(arg,lm.summ(f,v="Astar",lab="naive",conflevel=arg$conflevel))
}

validset.ols <- function(D,arg){
  f <- tryCatch.W.E(lm(paste0('Y~A+',paste(paste0('L.',1:arg$P),collapse='+')),data=D))
  data.frame(arg,lm.summ(f,v="A",lab="validset",conflevel=arg$conflevel))
}

mice.default.ols <- function(D,arg,nimp=10){
  D <- subset(D,select=-c(Ystd))
  imp <- tryCatch.W.E(mice(D,m=nimp, print = FALSE))
  fs <- with(imp$value,lm(as.formula(paste0('Y~A+',paste(paste0('L.',1:arg$P),collapse='+')))))
  data.frame(arg,mice.summ(fs,v="A",err=imp$warning,lab="impall",conflevel=arg$conflevel))
}

mice.noA.ols <- function(D,arg,nimp=10){
  D.noA <- subset(D,select=-c(A, Ystd))
  imp <- tryCatch.W.E(mice(D.noA,m=nimp, print = FALSE))
  fs <- with(imp$value,lm(as.formula(paste0('Y~Astar+',paste(paste0('L.',1:arg$P),collapse='+')))))
  data.frame(arg,mice.summ(fs,v="Astar",err=imp$warning,lab="impnoA",conflevel=arg$conflevel))
}


fiml2 <- function(D,arg){
  model <- paste0('Astar~A\nY~A+',paste(paste0('L.',1:arg$P),collapse='+'),
	'\nA~',paste(paste0('L.',1:arg$P),collapse='+'),'\nL.1~',paste(paste0('L.',2:arg$P),collapse='+'))
  f <- tryCatch.W.E(sem(model,data=D,missing="ml"))
  data.frame(arg,lavaan.summ(f,v="A",lab="fiml2",conflevel=arg$conflevel))
}

cc.mecor.ols <- function(D,arg){
  DcompL1 <- D[complete.cases(D[,"L.1"]),]
  f <- tryCatch.W.E(mecor(as.formula(paste0(' Y ~ MeasError(Astar, A) +',paste(paste0('L.',1:arg$P),collapse='+'))),, method = "rc_pooled1", data =   DcompL1, B = 0))
  data.frame(arg,mecor.summ(f,v="A",lab="ccmecor",conflevel=arg$conflevel))
}

imp.mecor.ols <- function(D,arg,nimp=10){
  D <- subset(D,select=-c(Ystd))
  Amis <- is.na(D$A)
  imp <- tryCatch.W.E(mice(D,m=nimp, print = FALSE))
  Dcomp <- complete(imp$value,"long")
  Dcomp$A[Dcomp$.id%in%which(Amis)] <- NA
 
  f <- list()
  OUT <- matrix(,ncol=2,nrow=nimp,dimnames=list(c(),c("est","estse"))) 
  for(im in 1:nimp){ 
    f[[im]] <- tryCatch.W.E(mecor(as.formula(paste0(' Y ~ MeasError(Astar, A) +',paste(paste0('L.',1:arg$P),collapse='+'))), method = "rc_pooled1", data =   Dcomp[Dcomp$.imp==im,], B = 0))
    OUT[im,"est"] <- as.numeric(f[[im]]$value$corfit$coefficients["A"])
    OUT[im,"estse"] <- sqrt(as.numeric(f[[im]]$value$corvar$var["A"]))
  }
  data.frame(arg,imp.summ(OUT,f,imp,lab="MImecor",conflevel=arg$conflevel))
}

bayes.ols <- function(D,arg,maxit=1000,burnin=500,thin=1,nchains=3){
  D <- subset(D,select=-c(Ystd))
  L <- grep('^L\\.[:0-9:]+$',colnames(D))
  data <- c(as.list(D[,-L]),L=list(D[,L,drop=FALSE]))
  data$n <- length(data$Y)
  data$P <- ncol(data$L)
  data$P1 <- 4L+data$P
  data$P2 <- data$P1+2L
  data$P3 <- data$P2+data$P-1L
  data$P4 <- data$P3+2L
  data$P5 <- data$P4+data$P-2L
  init <- list(tau=rep(1,4),b=rep(0,5L+3L*data$P))
  hyper <- list(alpha=1e-3,beta=1e-3,mu=0,nu=1e-4)
  modelstring <- paste0("model{
  for(i in 1:n){
  Astar[i]~dnorm(b[1]+b[2]*A[i],tau[1])
  Y[i]~dnorm(b[3]+b[4]*A[i]+inprod(b[5:P1],L[i,1:P]),tau[2])
  A[i]~dnorm(b[P2-1]+inprod(b[P2:P3],L[i,1:P]),tau[3])		
  L[i,1]~dnorm(b[P4-1]+inprod(b[P4:P5],L[i,2:P]),tau[4])
  }
  for(j in 1:4){
  tau[j]~dgamma(alpha,beta)
  }
  for(j in 1:",5L+3L*data$P,"){
  b[j]~dnorm(mu,nu)
  }
  }")
  model <- jags.model(textConnection(modelstring),data=c(data,hyper),inits=init,n.chains=nchains,quiet=T)
  update(model,n.iter=burnin,progress.bar="none")
  f <-tryCatch.W.E(coda.samples(model=model,variable.names=names(init),n.iter=maxit,thin=thin,progress.bar="none"))
  data.frame(arg,jags.summ(f,v="b[4]",lab="Bayes",conflevel=arg$conflevel,nchains=nchains))
}


################## execute all models ############################# 



## TEST ##
# setwd("C:/Users/mvansmeden/Dropbox/Projects-current/HG1_correcting_ME-MD-Conf/R")
# source("dependencies.R")
# DATA <- read.csv("results3/OUT.csv")
# arg <- DATA[10,]
# arg$conflevel=.95
#  D <- simdata(arg)
# 
#  system.time(arg <- naive.ols(D,arg))
#  system.time(arg <- validset.ols(D,arg))
#  system.time(arg <- mice.default.ols(D,arg))
#  system.time(arg <- fiml.default.ols(D,arg))
#  system.time(arg <- cc.mecor.ols(D,arg))
#  system.time(arg <- imp.mecor.ols(D,arg))
#  system.time(arg <- imp.mecor.ols(D,arg))
#  system.time(arg <- bayes.ols(D,arg,maxit=1e3,burnin=100,thin=1,nchains=3))
#  
