pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork8 && mpirun -np 8 zwork8
popd
