pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork10 && mpirun -np 10 zwork10
popd
