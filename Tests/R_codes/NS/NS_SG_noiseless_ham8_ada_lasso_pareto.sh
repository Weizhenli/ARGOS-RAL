#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1               # number of threads per rank
#SBATCH -n 25
#SBATCH --mem=220G 
#SBATCH -t 25:00:00

# try to control automatic multithreading
# ompthreads=1
# export OMP_NUM_THREADS=$ompthreads
# echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

export crit=AIC

# Commands to execute start here

cd /nobackup/cfzh32/for_r/SG/paper_test
# module load python
# module load r

R CMD BATCH NS_SG_noiseless_ham8_ada_lasso_pareto.R 
