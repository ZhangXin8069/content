pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork9 && mpirun -np 9 zwork9
popd
