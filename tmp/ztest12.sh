pushd /home/zhangxin/content/test
mpic++ ztest.cc -o ztest12 && mpirun -np 12 ztest12
popd
