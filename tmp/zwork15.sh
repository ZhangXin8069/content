pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork15 && mpirun.openmpi -np 15 zwork15
popd
