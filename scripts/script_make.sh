# source
_PATH=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
echo "PATH:"$_PATH
pushd ${_PATH}/../
source ./env.sh
popd

# init
_NAME=$(basename "$0")
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


# init
_NAME=$(basename "$0")
name='zwork'
work_name="bin"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# do
pushd ${tmp_path}
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
for i in $(seq 16); do
       echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${tmp_path}" >${name}${i}.sh
       echo "mpic++.openmpi ${work_path}/${name}.cc -o ${name}${i} && mpirun.openmpi -np ${i} ${name}${i}" >>${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd


# init
_NAME=$(basename "$0")
name='zexample'
work_name="bin"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# do
pushd ${tmp_path}
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
for i in $(seq 16); do
       echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${tmp_path}" >${name}${i}.sh
       echo "mpic++.openmpi ${work_path}/${name}.cc -o ${name}${i} && mpirun.openmpi -np ${i} ${name}${i}" >>${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd


# init
_NAME=$(basename "$0")
name='ztest'
work_name="bin"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# do
pushd ${tmp_path}
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
for i in $(seq 16); do
       echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${tmp_path}" >${name}${i}.sh
       echo "mpic++.openmpi ${work_path}/${name}.cc -o ${name}${i} && mpirun.openmpi -np ${i} ${name}${i}" >>${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd


# init
_NAME=$(basename "$0")
name='zrun'
work_name="bin"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# do
pushd ${tmp_path}
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
for i in $(seq 16); do
       echo "making ${name}${i}.sh in ${tmp_path}"
       echo "pushd ${tmp_path}" >${name}${i}.sh
       echo "mpic++.openmpi ${work_path}/${name}.cc -o ${name}${i} && mpirun.openmpi -np ${i} ${name}${i}" >>${name}${i}.sh
       echo "popd" >>${name}${i}.sh
done
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd
