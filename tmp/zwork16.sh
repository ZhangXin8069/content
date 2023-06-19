pushd /home/zhangxin/content/test
mpic++ zwork.cc -o zwork16 && mpirun -np 16 zwork16
popd
