pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample11 && mpirun.openmpi -np 11 zexample11
popd
