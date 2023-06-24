pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork10 && mpirun.openmpi -np 10 zwork10
popd
