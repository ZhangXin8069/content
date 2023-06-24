pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample4 && mpirun.openmpi -np 4 zexample4
popd
