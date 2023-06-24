pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample14 && mpirun.openmpi -np 14 zexample14
popd
