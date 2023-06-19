pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample2 && mpirun -np 2 zexample2
popd
