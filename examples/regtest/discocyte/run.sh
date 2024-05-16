#!/bin/bash 
#SBATCH --job-name=helf_2.50
#SBATCH --output=info_pip.log 
# 
#SBATCH -n 36
#SBATCH -N 1
#SBATCH --tasks-per-node=36
#SBATCH --time=12:00:00 
#SBATCH -p express3
#SBATCH --mail-type=ALL 
# create scr dir n copy input files
# Start the run
date
source /usr/local/gromacs/bin/GMXRC.bash
./nagata init.obj CHF
