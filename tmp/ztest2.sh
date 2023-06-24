pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest2 && mpirun.openmpi -np 2 ztest2
popd
