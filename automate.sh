#!/bin/bash
cd /home/anhvt89/Documents/lammps-28Jun14/bench-src-total-uncertainty-noabs/16x16x16-0.05
mpirun -np 4 ../lmp_openmpi < in.tensile.atomic_rad.txt;

cd /home/anhvt89/Documents/lammps-28Jun14/bench-src-ublb/16x16x16-0.05;
mpirun -np 4 ../lmp_openmpi < in.tensile.atomic_rad.txt;

# for i in /home/anhvt89/Documents/lammps-28Jun14/bench-src-total-uncertainty-abs/*/
# do
# 	echo ${i}
# 	cd ${i}
# 	mpirun -np 4 ../lmp_openmpi < in.tensile.atomic_rad.txt -log log.mpi.debug.tensile;
# 	echo "done ${i}"
# 	cd ..
# done

cd /home/anhvt89/Documents/lammps-28Jun14;