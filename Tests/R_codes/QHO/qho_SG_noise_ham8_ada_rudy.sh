#!/bin/bash

#SBATCH -p shared
#SBATCH -c 3               # number of threads per rank
#SBATCH -n 20
#SBATCH --mem=220G 
#SBATCH -t 72:00:00

# try to control automatic multithreading
# ompthreads=1
# export OMP_NUM_THREADS=$ompthreads
# echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Commands to execute start here
export OMP_NUM_THREADS=3
cd /nobackup/cfzh32/for_r/SG/paper_test
# module load python
# module load r

R CMD BATCH qho_SG_noise_ham8_ada_rudy.R
