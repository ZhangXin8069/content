pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork3 && mpirun.openmpi -np 3 zwork3
popd
