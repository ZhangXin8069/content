#!/bin/bash

#SBATCH --job-name=2pt_zx
#SBATCH --partition=64c512g
#SBATCH --output=test_2pt.sh.out
#SBATCH --error=test_2pt.sh.err
#SBATCH --array=0
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00

ulimit -s unlimited
ulimit -l unlimited

export OMP_NUM_THREADS=1
module load gcc/8.5.0
module load openmpi/4.1.1-gcc-8.5.0

list_id=list_all
prefix=$SLURM_SUBMIT_DIR

mkdir -p ini_out_file
mkdir -p log_file
mkdir -p Data

var=$(cat ${list_id})
for conf in ${var[*]}
do 
	./2pt.pl  ${conf} $prefix > ./ini_out_file/ini_${conf}.xml

	mpirun -mca btl_tcp_if_include ib0 -n 4  ../source_codes/chroma  -i ./ini_out_file/ini_${conf}.xml   -o ./ini_out_file/out_${conf}.xml >./log_file/log_${conf}   2>&1

	#  break
done


