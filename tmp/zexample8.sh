pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample8 && mpirun -np 8 zexample8
popd
