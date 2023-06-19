# source
_PATH=$(
    cd "$(dirname "$0")"
    pwd
)
echo "PATH:"$_PATH
pushd ${_PATH}/../
source ./env.sh
popd

# init
_NAME=$(basename "$0")
work_name="bin"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# do
pushd ${tmp_path}
rm scripts.sh
echo "###${_NAME} is running...:$(date "+%Y-%m-%d-%H-%M-%S")###"
echo "# >>> alias:$(date "+%Y-%m-%d-%H-%M-%S") >>>" >scripts.sh
for i in $(find ${tmp_path} -type f -name "*.sh"); do
    i=$(basename ${i})
    if [ ${i} = "scripts.sh" ];
    then
    continue
    fi
    echo "alias ${i}='bash ${tmp_path}/${i}' >>scripts.sh"
    echo "alias ${i}='bash ${tmp_path}/${i}'" >>scripts.sh
done
for i in $(find ${work_path} -type f -name "*.sh"); do
    i=$(basename ${i})
    if [ ${i} = "scripts.sh" ];
    then
    continue
    fi
    echo "alias ${i}='bash ${work_path}/${i}' >>scripts.sh"
    echo "alias ${i}='bash ${work_path}/${i}'" >>scripts.sh
done
echo "# <<< alias:$(date "+%Y-%m-%d-%H-%M-%S") <<<" >>scripts.sh
echo "###${_NAME} is done......:$(date "+%Y-%m-%d-%H-%M-%S")###"
popd

# done
