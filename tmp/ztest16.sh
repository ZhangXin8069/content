pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest16 && mpirun -np 16 ztest16
popd
