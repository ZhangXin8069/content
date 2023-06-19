pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest11 && mpirun -np 11 ztest11
popd
