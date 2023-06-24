pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork9 && mpirun.openmpi -np 9 zwork9
popd
