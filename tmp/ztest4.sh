pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest4 && mpirun -np 4 ztest4
popd
