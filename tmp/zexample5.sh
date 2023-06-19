pushd /home/zhangxin/content/test
mpic++ zexample.cc -o zexample5 && mpirun -np 5 zexample5
popd
