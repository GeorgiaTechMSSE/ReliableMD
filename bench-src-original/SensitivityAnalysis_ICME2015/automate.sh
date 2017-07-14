#!/bin/bash

cd /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_ICME2015/;

for i in /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_ICME2015/*/
do
	# go to local directory and cat the directory
	cd ${i}
	dir=$(basename $(pwd))

	# copy the original input
	cp /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_ICME2015/in.tensile.txt .

	# debug: 
	# echo "i = ${i}"
	# echo "dir = ${dir}"
	
	# change line 20
	sed -i.bak "20s|^.*$|pair_coeff	* * Al99.${dir}.eam.alloy Al|" in.tensile.txt
	
	# change line 63
	line63='fix def1 all print 20 "${p1} ${p2} ${p3} ${p4}" file '
	line63="$line63 Al99.${dir}.eam.alloy.def1.txt screen no" 
	sed -i.bak "63s|^.*$|${line63}|" in.tensile.txt

	# run the command
	mpirun -np 4 ../../lmp_openmpi < in.tensile.txt -log log.mpi.debug.tensile;
	echo "done ${i}"

	# copy the output upstair 
	cp Al99.${dir}.eam.alloy.def1.txt ..
	cd ..
done
