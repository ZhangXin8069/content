pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample14 && mpirun -np 14 zexample14
popd
