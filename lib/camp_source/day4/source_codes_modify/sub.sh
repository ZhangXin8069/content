#!/bin/bash

#SBATCH --job-name=make_source
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --exclusive
#SBATCH --output=make_source_code.out
#SBATCH --error=make_source_code.err
 
bash make.sh
