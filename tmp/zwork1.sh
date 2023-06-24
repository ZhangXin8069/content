pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork1 && mpirun.openmpi -np 1 zwork1
popd
