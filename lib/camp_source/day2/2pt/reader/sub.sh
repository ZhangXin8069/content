#!/bin/bash

#SBATCH --job-name=make_source
#SBATCH --partition=debug64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=task.out
#SBATCH --error=task.err

export PATH="/dssg/home/acct-phyww/phyww/qazhang/packages/anaconda3/bin:$PATH"
python read_iog.py > task.log
