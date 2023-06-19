pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork4 && mpirun -np 4 zwork4
popd
