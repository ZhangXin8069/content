pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork7 && mpirun.openmpi -np 7 zwork7
popd
