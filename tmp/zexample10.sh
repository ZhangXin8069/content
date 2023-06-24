pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample10 && mpirun.openmpi -np 10 zexample10
popd
