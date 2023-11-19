#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1               # number of threads per rank
#SBATCH -n 50
#SBATCH --mem=150G 
#SBATCH -t 10:00:00

# Commands to execute start here

cd /mnt/d/GitHub/ARGOS-RAL/Tests/R_codes/Keller
# module load python
# module load r

R CMD BATCH keller_segel_FD_noiseless_ham8_all_success_rate.R
