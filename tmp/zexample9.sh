pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample9 && mpirun.openmpi -np 9 zexample9
popd
