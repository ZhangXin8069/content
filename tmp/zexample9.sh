pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample9 && mpirun -np 9 zexample9
popd
