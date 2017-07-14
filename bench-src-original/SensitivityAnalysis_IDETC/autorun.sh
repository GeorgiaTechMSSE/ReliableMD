#!/bin/bash


for fid in `ls Al99*alloy`
do
	potentialFile=$(basename ${fid})
	# echo "${potentialFile} "

	catFRho=$( awk -F. -v OFS='\n' '{print $2}' <<< "$potentialFile")
	catRhoR=$( awk -F. -v OFS='\n' '{print $3}' <<< "$potentialFile")
	catPhiR=$( awk -F. -v OFS='\n' '{print $4}' <<< "$potentialFile")

	folderName=$(echo "${catFRho}.${catRhoR}.${catPhiR}")
	mkdir $folderName # create folder
	cd $folderName
	echo "currently in $(pwd)"

	cp /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_IDETC/in.tensile.txt .

	cp /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_IDETC/$potentialFile .

	# change line 20 - potentialFile input
	sed -i.bak "20s|^.*$|pair_coeff	* * ${potentialFile} Al|" in.tensile.txt
	
	# change line 63 - log output
	line63='fix def1 all print 20 "${p1} ${p2} ${p3} ${p4}" file ' # first part of line 63
	line63="$line63 ${potentialFile}.def1.txt screen no" # add to line 63
	sed -i.bak "63s|^.*$|${line63}|" in.tensile.txt

	# run the command
	mpirun -np 4 ../../lmp_openmpi < in.tensile.txt -log log.mpi.debug.tensile.txt;
	echo "done ${potentialFile} folder!"

	# copy result log file to database
	cp ${potentialFile}.def1.txt /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_IDETC/resultsDatabase/

	cd /home/anhvt89/Documents/lammps-28Jun14/bench-src-original/SensitivityAnalysis_IDETC; # go back to parent dir	
	# echo "${potentialFile} -> ${catFRho} -> ${catRhoR} -> ${catPhiR} " # debug
done
