pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample15 && mpirun.openmpi -np 15 zexample15
popd
