pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest8 && mpirun.openmpi -np 8 ztest8
popd
