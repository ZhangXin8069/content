pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork11 && mpirun.openmpi -np 11 zwork11
popd
