pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample13 && mpirun -np 13 zexample13
popd
