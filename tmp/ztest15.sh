pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest15 && mpirun.openmpi -np 15 ztest15
popd
