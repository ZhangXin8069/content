pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork3 && mpirun -np 3 zwork3
popd
