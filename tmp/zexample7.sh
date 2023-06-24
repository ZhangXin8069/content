pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample7 && mpirun.openmpi -np 7 zexample7
popd
