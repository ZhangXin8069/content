pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest6 && mpirun.openmpi -np 6 ztest6
popd
