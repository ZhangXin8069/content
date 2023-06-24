pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest9 && mpirun.openmpi -np 9 ztest9
popd
