pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest10 && mpirun -np 10 ztest10
popd
