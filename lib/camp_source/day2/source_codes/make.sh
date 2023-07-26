#!/bin/bash -l

INSTALL=/dssg/home/acct-phyww/phyww/chroma/chroma_double_gcc-8.5
export PATH=$INSTALL/bin:$PATH

module load openmpi/4.1.1-gcc-8.5.0
make clean
make -j 16

