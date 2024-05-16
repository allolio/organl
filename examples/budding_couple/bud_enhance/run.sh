#!/bin/bash 
#SBATCH --job-name=bud_enhance
#SBATCH --output=info_pip.log 
# 
#SBATCH -n 24 
#SBATCH -N 1
#SBATCH --tasks-per-node=24
#SBATCH --time=100:00:00 
#SBATCH -p nvidia
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<add_your_email_here>
#
date
source /usr/local/gromacs/bin/GMXRC.bash
./organl init.obj CHF LIP
date