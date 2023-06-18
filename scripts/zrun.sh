#!/bin/bash
#SBATCH -J test
#SBATCH --partition=cpueicc
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# mpirun -np 8 ~/a.out >result.log
mpirun -n 8 ~/a.out >result.log
