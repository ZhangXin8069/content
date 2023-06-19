pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork15 && mpirun -np 15 zwork15
popd
