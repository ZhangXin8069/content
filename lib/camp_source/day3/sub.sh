#!/bin/bash

#SBATCH --job-name=2pt_zx
#SBATCH --partition=64c512g
#SBATCH --output=test_2pt.sh.out
#SBATCH --error=test_2pt.sh.err
#SBATCH --array=0
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --exclusive

ulimit -s unlimited
ulimit -l unlimited

export OMP_NUM_THREADS=1
module load gcc/8.5.0
module load openmpi/4.1.1-gcc-8.5.0


	mpirun -mca btl_tcp_if_include ib0 -n 64  ./source_codes_Props/chroma  -geom 1 4 4 4  -i ./ini_4050.xml   -o ./out_4050.xml >./log_4050   2>&1



