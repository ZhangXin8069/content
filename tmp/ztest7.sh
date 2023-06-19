pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest7 && mpirun -np 7 ztest7
popd
