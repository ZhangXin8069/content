pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample5 && mpirun.openmpi -np 5 zexample5
popd
