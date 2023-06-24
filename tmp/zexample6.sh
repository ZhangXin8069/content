pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample6 && mpirun.openmpi -np 6 zexample6
popd
