pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest3 && mpirun -np 3 ztest3
popd
