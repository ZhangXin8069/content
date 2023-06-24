pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork4 && mpirun.openmpi -np 4 zwork4
popd
