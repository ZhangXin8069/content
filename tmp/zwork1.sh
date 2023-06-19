pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork1 && mpirun -np 1 zwork1
popd
