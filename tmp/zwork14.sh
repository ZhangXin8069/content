pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork14 && mpirun.openmpi -np 14 zwork14
popd
