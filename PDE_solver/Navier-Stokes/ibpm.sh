#!/bin/bash

# Request resources:
#SBATCH -p shared
#SBATCH -c 1 
#SBATCH -n 1          
#SBATCH --mem=5G      
#SBATCH --time=8:00:00



# Commands to be run:
cd /mnt/d/GitHub/ARGOS-RAL/pde_solver_data/Solve_NS/ibmp
./build/ibpm -ic von_karman/ibpm15000.bin -Re 100 -nx 450 -ny 200 -ngrid 4 -length 9 -xoffset -1 -yoffset -2 -xshift 0.75 -nsteps 1510 -restart 100 -tecplot 10 -outdir try/ -geom von_karman/cylinder.geom -dt 0.02 -ubf 0 -tecplotallgrids 0