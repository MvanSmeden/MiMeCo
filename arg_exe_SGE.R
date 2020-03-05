args <- commandArgs(TRUE)
i <- as.integer(args[1L])
B <- as.integer(args[2L])

setwd('') # set work directory
filepath <- '' # set file path

source("data_and_parameter_functions.R")
source("data_gen_and_check.R")
source("models.R")

dr <- paste0(filepath,i,'.csv')
gendgm_and_check((i-1)*B+1,dr,append=FALSE)
for(j in 2:B) gendgm_and_check((i-1)*B+j,dr,append=TRUE)