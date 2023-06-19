pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample4 && mpirun -np 4 zexample4
popd
