pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest5 && mpirun.openmpi -np 5 ztest5
popd
