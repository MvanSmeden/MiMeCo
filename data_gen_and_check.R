################################################################
############ GENERATE DATA ANC CHECKS ##########################
################################################################

###############################################
# Author: M van Smeden                        #
# Date first version: March 15 2019           #
# Simulations confounding, measurement error  #
# and missing values                          #
###############################################

################## generate new  ##################################
gendgm_and_check <-function(iter,filepath,append=TRUE){
  seed1 <- as.numeric(Sys.time())+iter*10
  set.seed(seed1)
    arg <- sampl_sim_cond(bd_N=c(100,2000),bd_P=c(2,6),bd_beta=10,bd_R2=c(0.10,0.30),bd_rho=c(0.10,0.40),
                        bd_delta=c(0.50,0.99),bd_kappa=c(0.50,0.99),bd_lambda=c(0.50,0.99),
                        bd_validfrac = c(0.10,0.30),bd_prop.mis=0.30) # RANGES OF SIMULATION FACTORS
      arg$seed.dgm <- seed1
        dgmtime <- system.time(arg<-dgmpar_sim_cond(arg))
          arg$dgmtime_sec <- as.numeric(dgmtime["elapsed"])
        lsctime <- system.time(arg<-large_sample_check(arg))
          arg$lsctime_sec <- as.numeric(lsctime["elapsed"])
            arg$timing <- Sys.time()

      data.table::fwrite(arg,file=filepath,append=append)
  cat(paste0('\r',iter, '; this is row number ',length(count.fields(filepath)))); flush.console()
}

################### function large sample check ###################

large_sample_check <- function(arg,lsc_n=1e5){
  arg_val <- arg
  arg_val$N <- lsc_n
  DATA <- simdgm(arg_val)
  fit.compl <- lm(Y~.,data=DATA)
  arg$lsc.beta <- as.numeric(fit.compl$coefficients["A"])
  arg$lsc.R2 <- as.numeric(summary(fit.compl)$r.squared)
  fit.conf <- lm(Y~A,data=DATA)
  arg$lsc.delta <- as.numeric(fit.conf$coefficients["A"])/arg$beta

  DATA$Astar <- DATA$A+rnorm(arg_val$N,mean=0,sd=arg$tau)
  fit.me <- lm(as.formula(paste0('Y~Astar+',paste(paste0('L.',1:arg$P),collapse='+'))),data=DATA)
  arg$lsc.lambda <- as.numeric(fit.me$coefficients["Astar"])/arg$beta
  arg$lsc.cor.A.Astar <- cor(DATA$Astar,DATA$A)

  DATA$Ystd <-DATA[,"Y"]/sqrt(arg$var_Y) 
  dmd <- model.matrix(as.formula(arg$mddesmat),DATA)
  mdp <- plogis(dmd%*%matrix(c(arg$phi,rep(arg$omega,ncol(dmd)-1)),ncol=1))  
  arg$lsc.prop.mis <- mean(mdp)
  R <- 1*(runif(arg_val$N)< mdp)
  fit.md <- lm(as.formula(paste0('Y~A+',paste(paste0('L.',1:arg$P),collapse='+'))),data=DATA[R==0,])
  arg$lsc.kappa <- as.numeric(fit.md$coefficients["A"])/arg$beta
  
  fit.me.md <- lm(as.formula(paste0('Y~Astar+',paste(paste0('L.',1:arg$P),collapse='+'))),data=DATA[R==0,])
  arg$lsc.beta.me.md <- as.numeric(fit.me.md$coefficients["Astar"])
  fit.me.md.conf <- lm(Y~Astar,data=DATA[R==0,])
  arg$lsc.beta.me.md.conf <- as.numeric(fit.me.md.conf$coefficients["Astar"])

  arg
}


############### function data sample check ##################
sim_sample_check <- function(D,arg){
  arg$D.Nmis <-  sum(is.na(D$L.1))
  arg$D.Nvalid <-  sum(!is.na(D$L.1))
  arg$D.val.A <- sum(!is.na(D$A))  
  arg$D.var.Y <- var(D$Y)
  arg
}

#### END ###