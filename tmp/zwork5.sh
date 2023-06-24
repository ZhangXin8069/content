pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork5 && mpirun.openmpi -np 5 zwork5
popd
