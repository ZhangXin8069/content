pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork12 && mpirun.openmpi -np 12 zwork12
popd
