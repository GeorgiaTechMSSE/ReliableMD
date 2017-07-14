#!/bin/bash
# quick run after compile

reset;
make openmpi;

mv lmp_openmpi ../bench/lmp_openmpi_cxsc;
cd ../bench;
# # rm -rf dump.*.txt;
mpirun -v -np 4 lmp_openmpi_cxsc < in.tensile.atomic_rad.txt;
# mv dump*.txt dump/;
# # diff log.debug.tensile log.ref.tensile;
cd ../src-cxsc;