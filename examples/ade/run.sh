#!/bin/bash 
#SBATCH --job-name=ade
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
./organl init.obj ADE
date