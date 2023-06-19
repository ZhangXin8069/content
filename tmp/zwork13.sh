pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork13 && mpirun -np 13 zwork13
popd
