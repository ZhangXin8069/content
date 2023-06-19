pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork14 && mpirun -np 14 zwork14
popd
