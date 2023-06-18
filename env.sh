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
source ${_HOME}/include/.bashrc

# do
## export
export CPATH=$CPATH:$HOME/lib/
export CPATH=$CPATH:$HOME/lib/eigen-3.4.0/
export CPATH=$CPATH:$HOME/lib/openmpi-4.1.2/
export PATH=$PATH:$HOME/cling/bin/
# export PYTHONPATH=/home/aistudio/external-libraries:$PYTHONPATH
# export LD_LIBRARY_PATH=/home/aistudio/external-libraries/quda/build/lib/libquda.so:$LD_LIBRARY_PATH

## add alias
alias noita="pushd /home/zhangxin/Game/Noita\ v20230311 && wine noita.exe && popd"
alias dwarf="pushd /home/zhangxin/Game/Dwarf\ Fortress && wine Dwarf\ Fortress.exe && popd"
alias rain="pushd /home/zhangxin/Game/Rain\ World\ v1.9.07b && wine RainWorld.exe && popd"
alias oriwotw="pushd /home/zhangxin/Game/Ori\ and\ the\ Will\ of\ the\ Wisps && wine oriwotw.exe && popd"
alias deadcells="pushd /home/zhangxin/Game/Dead\ Cells && wine deadcells.exe  && popd"
alias space="pushd /home/zhangxin/Packages && wine SpaceSniffer.exe  && popd"
alias winrar="pushd /home/zhangxin/Packages/WinRARPortable && wine WinRARPortable.exe && popd"
alias matlab="pushd /home/zhangxin/MATLAB/R2023b/bin && bash matlab && popd"

# done