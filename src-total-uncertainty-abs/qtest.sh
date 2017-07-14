#!/bin/bash
reset;
cp pair.cpp.v2 pair.cpp
cp pair.h.v2 pair.h
cp pair_eam_rad.cpp.v2 pair_eam_rad.cpp
cp pair_eam_rad.h.v2 pair_eam_rad.h
cp verlet.cpp.v2 verlet.cpp
cp atom.cpp.v2 atom.cpp
cp atom.h.v2 atom.h
cp atom_vec_atomic_rad.cpp.v2 atom_vec_atomic_rad.cpp
cp atom_vec_atomic_rad.h.v2 atom_vec_atomic_rad.h
make openmpi;