#!/bin/bash
reset;

# copy and compile new file here
# cp atom.cpp.v3 atom.cpp
# cp atom.h.v3 atom.h
# cp atom_vec_atomic_rad.cpp.v3 atom_vec_atomic_rad.cpp
# cp atom_vec_atomic_rad.h.v3 atom_vec_atomic_rad.h
# cp compute.h.v3.3 compute.h
# cp compute_pressure.cpp.v3.3 compute_pressure.cpp
# cp compute_pressure.h.v3.3 compute_pressure.h
# cp compute_temp.cpp.v3.3 compute_temp.cpp
# cp compute_temp.h.v3.3 compute_temp.h
# cp fix_nh.cpp.v3 fix_nh.cpp
# cp fix_nh.h.v3 fix_nh.h
# cp pair.cpp.v3 pair.cpp
# cp pair.h.v3 pair.h
cp pair_eam_rad.cpp.v3 pair_eam_rad.cpp
# cp pair_eam_rad.h.v3 pair_eam_rad.h
# cp thermo.cpp.v3 thermo.cpp
# cp thermo.h.v3 thermo.h
# cp verlet.cpp.v3 verlet.cpp
make openmpi

mv lmp_openmpi ../bench/lmp_openmpi_mod;
cd ../bench;
# rm -rf dump.*.txt;
mpirun -v -np 4 lmp_openmpi_mod < in.tensile.atomic_rad.txt;
# mv dump*.txt dump/;
# diff log.debug.tensile log.ref.tensile;
cd ../src;