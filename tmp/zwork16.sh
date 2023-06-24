pushd /home/zhangxin/content/test
mpic++.openmpi zwork.cc -o zwork16 && mpirun.openmpi -np 16 zwork16
popd
