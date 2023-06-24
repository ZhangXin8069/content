pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample1 && mpirun.openmpi -np 1 zexample1
popd
