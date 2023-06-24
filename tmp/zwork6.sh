pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork6 && mpirun.openmpi -np 6 zwork6
popd
