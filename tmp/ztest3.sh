pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest3 && mpirun.openmpi -np 3 ztest3
popd
