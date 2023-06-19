pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork5 && mpirun -np 5 zwork5
popd
