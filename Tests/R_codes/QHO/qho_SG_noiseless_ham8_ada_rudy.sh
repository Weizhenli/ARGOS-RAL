#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1               # number of threads per rank
#SBATCH -n 25
#SBATCH --mem=200G 
#SBATCH -t 40:00:00

# try to control automatic multithreading
# ompthreads=1
# export OMP_NUM_THREADS=$ompthreads
# echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS
# Commands to execute start here
export OMP_NUM_THREADS=1
cd /nobackup/cfzh32/for_r/SG/paper_test
# module load python
# module load r

R CMD BATCH qho_SG_noiseless_ham8_ada_rudy.R
