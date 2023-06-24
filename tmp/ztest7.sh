pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest7 && mpirun.openmpi -np 7 ztest7
popd
