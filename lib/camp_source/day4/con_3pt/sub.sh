#!/bin/bash

#SBATCH --job-name=3pt_hj
#SBATCH --partition=64c512g
#SBATCH --output=test_3pt.sh.out
#SBATCH --error=test_3pt.sh.err
#SBATCH --array=0
#SBATCH -N 1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=02:20:00
#SBATCH --exclusive

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
mkdir -p Jobs

var=$(cat ${list_id})

if [ ${SLURM_ARRAY_TASK_ID} -eq "0" ]; then
        for i in ${var[*]} #loop list
        do
              mkdir -p "./Jobs/${i}" #make job list for each conf, delete a job after running one conf
        done
else
        sleep $(5+(${SLURM_ARRAY_TASK_ID}*1)) #sleep time
fi


while true
do
        count=`ls ./Jobs|wc -w`
        if [ $count -eq "0"  ]; then
                break
        fi
        i=`ls -r ./Jobs | tail -1` #ls job list, pick last one
        rm -rf "./Jobs/${i}" #delet one job, after running one conf

	./3pt.pl  ${i} $prefix > ./ini_out_file/ini_${i}.xml
	HIP_DB=0 QUDA_ENABLE_P2P=0 QUDA_ENABLE_TUNING=0 \
	mpirun -mca btl_tcp_if_include ib0 -n 64  ./chroma  -i ./ini_out_file/ini_${i}.xml   -o ./ini_out_file/out_${conf}.xml >./log_file/log_${i}   2>&1
	#break
done


