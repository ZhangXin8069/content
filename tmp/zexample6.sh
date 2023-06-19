pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample6 && mpirun -np 6 zexample6
popd
