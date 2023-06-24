pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest14 && mpirun.openmpi -np 14 ztest14
popd
