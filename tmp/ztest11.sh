pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest11 && mpirun.openmpi -np 11 ztest11
popd
