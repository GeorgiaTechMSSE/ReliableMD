#!/bin/bash

# cp fix_nh.cpp.bak.21Nov15 fix_nh.cpp;
# cp fix_nh.h.bak.11Jan16 fix_nh.h;

# reset;
clear;
make openmpi;

mv lmp_openmpi ../bench/lmp_openmpi_cxsc;
cd ../bench;
# # rm -rf dump.*.txt;
mpirun -v -np 4 lmp_openmpi_cxsc < in.tensile.atomic_rad.txt;
# mv dump*.txt dump/;
# # diff log.debug.tensile log.ref.tensile;
cd ../src-cxsc;