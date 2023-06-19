pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample3 && mpirun -np 3 zexample3
popd
