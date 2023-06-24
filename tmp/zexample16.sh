pushd /home/zhangxin/content/test
mpic++.openmpi zexample.cc -o zexample16 && mpirun.openmpi -np 16 zexample16
popd
