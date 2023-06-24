pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample8 && mpirun.openmpi -np 8 zexample8
popd
