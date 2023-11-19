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

# Commands to execute start here

cd /mnt/d/GitHub/ARGOS-RAL/Tests/R_codes/Heat
# module load python
# module load r

R CMD BATCH heat_SG_noise_ham8_ada_lasso_pareto_success_rate.R 
