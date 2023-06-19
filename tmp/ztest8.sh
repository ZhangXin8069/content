pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest8 && mpirun -np 8 ztest8
popd
