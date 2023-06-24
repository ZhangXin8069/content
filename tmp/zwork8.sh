pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork8 && mpirun.openmpi -np 8 zwork8
popd
