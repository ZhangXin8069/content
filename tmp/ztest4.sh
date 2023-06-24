pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest4 && mpirun.openmpi -np 4 ztest4
popd
