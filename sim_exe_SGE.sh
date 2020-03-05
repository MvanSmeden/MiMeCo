#!/bin/bash
#$ -S /bin/bash
#$ -q all.q
#$ -N sim_MvS
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
R CMD BATCH --vanilla "--args $i $B" 
path/sim_exe_SGE.R 
path/sim_Rout_SGE/$i.Rout