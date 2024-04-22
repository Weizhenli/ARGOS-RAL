#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1             # number of threads per rank
#SBATCH -n 50
#SBATCH --mem=100G 
#SBATCH -t 05:00:00

# try to control automatic multithreading
ompthreads=1
export OMP_NUM_THREADS=$ompthreads
echo 'OMP_NUM_THREADS is ' $OMP_NUM_THREADS

# Commands to execute start here

cd /mnt/d/GitHub/ARGOS-RAL/Tests/R_codes/Noise # cd to the directory with the R codes
# module load python
# module load r

# time R CMD BATCH bur_SG_noise_ham8.R
time R CMD BATCH random_noise_test4.R random_noise_test_p${p}_d${d}_sd${sd}_opt.Rout
