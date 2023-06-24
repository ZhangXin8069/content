pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest12 && mpirun.openmpi -np 12 ztest12
popd
