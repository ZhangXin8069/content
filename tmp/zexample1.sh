pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample1 && mpirun -np 1 zexample1
popd
