pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample2 && mpirun.openmpi -np 2 zexample2
popd
