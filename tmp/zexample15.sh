pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample15 && mpirun -np 15 zexample15
popd
