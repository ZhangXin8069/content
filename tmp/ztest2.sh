pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest2 && mpirun -np 2 ztest2
popd
