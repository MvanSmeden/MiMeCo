#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N args_MvS
#$ -hard
#$ -cwd
#$ -e /exports/
#$ -o /exports/
#$ -j Y
#$ -V
#$ -m a
#$ -M email

##echo 'hello'
module load R/3.5.0
R CMD BATCH --vanilla "--args $i $B" 
path/arg_exe_SGE.R 
path/$i.Rout