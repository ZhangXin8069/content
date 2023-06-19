pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork11 && mpirun -np 11 zwork11
popd
