pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest15 && mpirun -np 15 ztest15
popd
