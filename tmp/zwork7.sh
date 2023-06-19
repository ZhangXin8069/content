pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork7 && mpirun -np 7 zwork7
popd
