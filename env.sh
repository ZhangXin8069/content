# init
_HOME=$(
    pwd
)
echo 'HOME:'${_HOME}
_NAME=$(basename "$0")
name='test'
work_name="test"
tmp_name="tmp"
work_path=${_HOME}/${work_name}
tmp_path=${_HOME}/${tmp_name}

# source
## mkdir
mkdir ${_HOME}/bin
mkdir ${_HOME}/include
mkdir ${_HOME}/lib
mkdir ${_HOME}/scripts
mkdir ${_HOME}/test
mkdir ${_HOME}/tmp

source ${_HOME}/tmp/scripts.sh

# do
## export
export CPATH=$CPATH:$HOME/lib/
export CPATH=$CPATH:$HOME/lib/eigen-3.4.0/
export CPATH=$CPATH:$HOME/lib/openmpi-4.1.2/
export PATH=$PATH:$HOME/cling/bin/
# export PYTHONPATH=/home/aistudio/external-libraries:$PYTHONPATH
# export LD_LIBRARY_PATH=/home/aistudio/external-libraries/quda/build/lib/libquda.so:$LD_LIBRARY_PATH

## add alias

# done
