pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest10 && mpirun.openmpi -np 10 ztest10
popd
