pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork13 && mpirun.openmpi -np 13 zwork13
popd
