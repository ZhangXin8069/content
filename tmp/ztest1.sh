pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest1 && mpirun.openmpi -np 1 ztest1
popd
