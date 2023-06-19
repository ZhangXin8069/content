# source
_PATH=`pushd "$(dirname "$0")";pwd`
pushd ${_PATH}/../
source ./env.sh
popd

# init
_NAME=`basename "$0"`
name='work'
work_name="test"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# do
pushd ${tmp_path}
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
for ((i = 0; i < 10; i++)); do
    echo "making ${name}${i}.sh in ${tmp_path}"
    echo "pushd ${tmp_path}" >${name}${i}.sh
    echo "nvcc -o ${name}${i} ${work_path}/${name}${i}.cu && ./${name}${i}" >>${name}${i}.sh
    echo "popd" >>${name}${i}.sh
done
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd

pushd ${tmp_path}
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
name="zwork"
for i in $(seq 16)
do
echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${work_path}" > ${name}${i}.sh
       echo "mpic++ ${name}.cc -o ${name}${i} && mpirun -np ${i} ${name}${i}" >> ${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
name="zexample"
for i in $(seq 16)
do
echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${work_path}" > ${name}${i}.sh
       echo "mpic++ ${name}.cc -o ${name}${i} && mpirun -np ${i} ${name}${i}" >> ${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
name="ztest"
for i in $(seq 16)
do
echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${work_path}" > ${name}${i}.sh
       echo "mpic++ ${name}.cc -o ${name}${i} && mpirun -np ${i} ${name}${i}" >> ${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd

# done
