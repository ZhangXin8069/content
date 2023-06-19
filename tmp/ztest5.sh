pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest5 && mpirun -np 5 ztest5
popd
