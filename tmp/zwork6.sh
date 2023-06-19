pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork6 && mpirun -np 6 zwork6
popd
