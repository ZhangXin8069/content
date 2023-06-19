pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest13 && mpirun -np 13 ztest13
popd
