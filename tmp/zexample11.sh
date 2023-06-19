pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample11 && mpirun -np 11 zexample11
popd
