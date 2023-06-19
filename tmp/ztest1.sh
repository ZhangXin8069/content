pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest1 && mpirun -np 1 ztest1
popd
