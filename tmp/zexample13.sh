pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample13 && mpirun.openmpi -np 13 zexample13
popd
