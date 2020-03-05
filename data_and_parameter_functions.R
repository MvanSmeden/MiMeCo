################################################################
############ SIMULATION CONDITIIONS FUNCTIONS ##################
################################################################

###############################################
# Author: M van Smeden                        #
# Date first version: March 15 2019           #
# Simulations confounding, measurement error  #
# and missing values                          #
###############################################

################## data generating functions ##################
simdgm <- function(arg){
  S <- matrix(arg$rho,ncol=arg$P+1,nrow=arg$P+1)
    diag(S) <- 1
    X <- MASS::mvrnorm(arg$N,mu=rep(0,arg$P+1),Sigma=S)
    f <- matrix(X[,1]*arg$beta,ncol=1)+X[,-1]%*%matrix(rep(arg$gamma,arg$P),ncol=1)
    Y <- f + rnorm(arg$N,mean=0,sd=arg$sigma)
  data.frame(Y=Y,A=X[,1],L=X[,-1])
}

simdata <- function(arg){
  D <- simdgm(arg)
  D$Astar <- D$A+rnorm(arg$N,mean=0,sd=arg$tau)
  D$Ystd <- D$Y/sqrt(arg$var_Y)
    Dmd <- model.matrix(as.formula(as.character(arg$mddesmat)),D)
     mdp <- plogis(Dmd%*%matrix(c(arg$phi,rep(arg$omega,ncol(Dmd)-1)),ncol=1))
        R <- (runif(arg$N) < mdp)
          D[R,"L.1"] <- NA
    D$Astar <- D$A+rnorm(arg$N,mean=0,sd=arg$tau)
    nV <-  sample(1:arg$N,size=arg$N-round(arg$validfrac*arg$N),replace=F)
      D[1:arg$N%in%nV,"A"] <- NA
  D
}

################## function  dgm ###############################
find_gamma_sigma <- function(params,arg){
  gamma <- params[1] 
  Sigma <- params[2]
    var_Ls  <- arg$P+arg$rho*(arg$P-1)*arg$P
    var_A <- 1
    cov_A_Ls <- arg$P*arg$rho
    theta <- cov_A_Ls/var_Ls
    cov_Y_Ls <- gamma*var_Ls + arg$beta*var_A*cov_A_Ls 
    cov_Y_A  <- arg$beta*var_A + theta*gamma*var_Ls 
    bias <- gamma*theta*(var_Ls/var_A)
    var_beta_A <- arg$beta^2
    var_gamma_Ls <- gamma^2*var_Ls
    var_f <- arg$beta^2 + gamma^2*arg$P + gamma^2*arg$rho*(arg$P-1)*arg$P+ 2*arg$beta*gamma*arg$P*arg$rho
    var_Y <- var_f + Sigma^2
    Rsq <- var_f/var_Y
    targ <- arg$beta*arg$delta
  log(abs(Rsq-arg$R2)+(abs(bias-(targ-arg$beta))))
}

################## function  measurement error ###################
find_wvar <- function(param,arg){
  tau <- param
    cov_A_Ls <- arg$P*arg$rho
    var_Ls <- arg$P+arg$rho*(arg$P-1)*arg$P
    var_A <- 1
    atten <- 1/(1+tau^2)
    Rsp <- cov_A_Ls^2/(var_Ls*(var_A+tau^2))
    attfac <- (atten-Rsp)/(1-Rsp)
  log(abs(attfac-arg$lambda))
}

################## function  missing data ########################
find_amp_dgm <- function(params,DATA,Dmd,U,arg){
  phi <- params[1] 
  omega <- params[2]
  mdp <- plogis(Dmd%*%matrix(c(phi,rep(omega,ncol(Dmd)-1)),ncol=1))  
  R <- 1*(U < mdp)
  est <- lm(paste0('Y~A+',paste(paste0('L.',1:arg$P),collapse='+')),data=DATA[R==0,])$coef['A']
  #est <- coef(.lm.fit(x=data.matrix(DATA[R==0,2:(arg$P+2)]),y=data.matrix(DATA[R==0,1])))[1]
  pR <- mean(mdp)
  log(abs(pR-arg$prop.mis) + abs(est/arg$beta-arg$kappa))
}

sim_amp_dgm <- function(arg,nsim=1e+05){
  argsim <- arg
  argsim$N <- nsim
  DATA <- simdgm(argsim)
  DATA$Astar <- DATA$A+rnorm(argsim$N,mean=0,sd=argsim$tau)
  DATA$Ystd <- DATA$Y/sqrt(argsim$var_Y)
  Dmd <- model.matrix(as.formula(as.character(argsim$mddesmat)),DATA)
  U <- runif(nsim)
  optim(par=c(-1,.2),fn=find_amp_dgm,DATA=DATA,Dmd=Dmd,U=U,arg=argsim, control = list(reltol = 1e-10))
}


################## function generate sim condition #############
sampl_sim_cond <- function(bd_N,bd_P,...){
  argg <- as.list(sys.call())[-1]
  nargg <- length(argg)
  TMP <- matrix(,nrow=1,ncol= nargg,dimnames=list(c(),substring(names(argg), 4)))
  for(j in 1:2){
    values <-  lazyeval::lazy_eval(toString(argg[j]))
    TMP[1,j] <- sample(min(values):max(values),1)
  }
  for(j in 3:nargg){
    values <-  lazyeval::lazy_eval(toString(argg[j]))
    TMP[1,j] <- runif(1,min(values),max(values))
  }
  data.frame(TMP)  
}

dgmpar_sim_cond <- function(arg){ 
  dgmparam <- optim(par=c(-10,100),fn=find_gamma_sigma,arg=arg)
  arg$gamma <- dgmparam$par[1]
  arg$sigma <- abs(dgmparam$par[2])
  arg$var_Y <- arg$beta^2+arg$gamma^2*arg$P+arg$gamma^2*arg$rho*(arg$P-1)*arg$P+2*arg$beta*arg$gamma*arg$P*arg$rho+arg$sigma^2
  
  meparam <- optimize(f=find_wvar,c(0,10^5),arg)
  arg$tau <- meparam$minimum

  arg$mddesmat <- "~Ystd * Astar + L.2"
  mdparam <- sim_amp_dgm(arg)
  arg$phi <-  mdparam$par[1]
  arg$omega <-  mdparam$par[2]
 
  arg$optval.dgm <- dgmparam$value
  arg$optval.me <- meparam$objective
  arg$optval.md <- mdparam$value
  
  arg
}

#### END ###