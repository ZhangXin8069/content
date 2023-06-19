pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest14 && mpirun -np 14 ztest14
popd
