pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample16 && mpirun -np 16 zexample16
popd
