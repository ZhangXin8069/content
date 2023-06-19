pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample12 && mpirun -np 12 zexample12
popd
