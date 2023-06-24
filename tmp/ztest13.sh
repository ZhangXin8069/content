pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest13 && mpirun.openmpi -np 13 ztest13
popd
