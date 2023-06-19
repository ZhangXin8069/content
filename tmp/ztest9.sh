pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest9 && mpirun -np 9 ztest9
popd
