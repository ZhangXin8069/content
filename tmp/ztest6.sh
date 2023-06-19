pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest6 && mpirun -np 6 ztest6
popd
