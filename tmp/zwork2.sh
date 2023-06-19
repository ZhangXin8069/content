pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork2 && mpirun -np 2 zwork2
popd
