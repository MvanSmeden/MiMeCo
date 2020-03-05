#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N proc_MvS
#$ -hard
#$ -cwd
#$ -e path/sim_Rout_SGE
#$ -o path/sim_Rout_SGE
#$ -j Y
#$ -V
#$ -m a
#$ -M email

##echo 'hello'
module load R/3.5.0
R CMD BATCH --vanilla "--args $B $L $U $X" 
path/process.R 
path/proc.Rout
