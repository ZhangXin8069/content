pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample12 && mpirun.openmpi -np 12 zexample12
popd
