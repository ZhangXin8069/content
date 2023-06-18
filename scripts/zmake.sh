path="~/content/zcc"
name="zwork"
for i in $(seq 16)
do
       echo "cd ${path}" > ${name}${i}.sh
       echo "mpic++ ${name}.cc -o ${name}${i} && mpirun -np ${i} ${name}${i}" >> ${name}${i}.sh
done
name="zexample"
for i in $(seq 16)
do
       echo "cd ${path}" > ${name}${i}.sh
       echo "mpic++ ${name}.cc -o ${name}${i} && mpirun -np ${i} ${name}${i}" >> ${name}${i}.sh
done
name="ztest"
for i in $(seq 16)
do
       echo "cd ${path}" > ${name}${i}.sh
       echo "mpic++ ${name}.cc -o ${name}${i} && mpirun -np ${i} ${name}${i}" >> ${name}${i}.sh
done
