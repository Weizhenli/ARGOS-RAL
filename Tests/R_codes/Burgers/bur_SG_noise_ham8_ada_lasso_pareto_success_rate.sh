#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1               # number of threads per rank
#SBATCH -n 12
#SBATCH --mem=60G 
#SBATCH -t 24:00:00

# try to control automatic multithreading
# ompthreads=1
# export OMP_NUM_THREADS=$ompthreads
# echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

export crit=AIC

# Commands to execute start here

cd /nobackup/cfzh32/for_r/SG/paper_test
# module load python
# module load r

R CMD BATCH bur_SG_noise_ham8_ada_lasso_pareto_success_rate.R 
