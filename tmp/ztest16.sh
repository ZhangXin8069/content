pushd /home/zhangxin/content/test
mpic++.openmpi ztest.cc -o ztest16 && mpirun.openmpi -np 16 ztest16
popd
