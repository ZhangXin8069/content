pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample3 && mpirun.openmpi -np 3 zexample3
popd
