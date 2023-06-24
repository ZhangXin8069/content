pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork2 && mpirun.openmpi -np 2 zwork2
popd
