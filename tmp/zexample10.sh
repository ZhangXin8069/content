pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample10 && mpirun -np 10 zexample10
popd
