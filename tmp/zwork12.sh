pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork12 && mpirun -np 12 zwork12
popd
