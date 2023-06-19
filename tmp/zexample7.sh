pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample7 && mpirun -np 7 zexample7
popd
