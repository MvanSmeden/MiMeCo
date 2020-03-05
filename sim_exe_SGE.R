args <- commandArgs(TRUE)
i <- as.integer(args[1L])
B <- as.integer(args[2L])
I <- (i-1)*B+1
I <- I+0:(B-1L)

setwd('path') # change path

.libPaths('path/packages_SGE') #change path
library(mice)
library(lavaan)
library(mecor)
library(rjags)
source("data_and_parameter_functions.R")
source("data_gen_and_check.R")
source("models.R")
load("/exports/clinicalepi/Bas/MvS/args.RData")

Models <- c("gs.ols","ccunivar.ols","univar.ols","naive.ols",
	"validset.ols","mice.default.ols","mice.noA.ols",
	"fiml2","imp.mecor.ols","cc.mecor.ols","bayes.ols"
)

filepath <- # add path

N_RES <- 168

single_run <- function(i,filename,append,models=Models){
	Sta = Sys.time()
	arg <- results[i,]
	set.seed(arg$seed.dgm)
	D <- simdata(arg)
	arg$conflevel <- 0.95
	arg <<- arg
	for(j in 1:length(models)){
		#arg <- do.call(models[j],args=list(D=D,arg=arg))
		tmp <- tryCatch.W.E(do.call(models[j],args=list(D=D,arg=arg)))
		if(!'error'%in%class(tmp$value)) arg <- tmp$value else{
			arg <- data.frame(arg,as.list(rep(NA,N_RES-ncol(arg))))
			break
		}
	}
	End = Sys.time()
	arg$modeltime_sec <-abs(round(as.numeric(difftime(time1=Sta,time2=End,units="secs")),3))
	data.table::fwrite(arg,file=filename,append=append)
}

single_run(I[1],paste0(filepath,i,'.csv'),append=FALSE)
for(k in 2:B) single_run(I[k],paste0(filepath,i,'.csv'),append=TRUE)