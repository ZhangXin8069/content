#!/bin/bash 

#SBATCH --job-name=conf_gen
#SBATCH --mail-user=wei.wang@sjtu.edu.cn
#SBATCH --mail-type=end
#SBATCH --partition=a100
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --exclude=gpu12
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
##SBATCH --exclusive
##SBATCH --time=24:00:00
#SBATCH --output=hmc_%j.out
#SBATCH --error=hmc_%j.err

ulimit -s unlimited
ulimit -l unlimited

export OMP_NUM_THREADS=8
module purge
module load gcc/8.3.1
module load openmpi/4.1.1-gcc-8.3.1
module load cuda/11.5.0   





#QUDA_ENABLE_P2P=3 QUDA_ENABLE_TUNING=1 QUDA_RESOURCE_PATH=./ \
    mpirun -n 1 ./hmc  -i  hmc.prec_wilson.ini.xml  -o mout -geom 1 1 1 1 
