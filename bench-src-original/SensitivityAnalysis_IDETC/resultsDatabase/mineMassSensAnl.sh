#!/bin/bash

minStressFile=minMassSensStress.def.txt;
maxStressFile=maxMassSensStress.def.txt;
numberOfLines=$(wc -l Al99.NegFRho.ZeroRhoR.ZeroPhiR.eam.alloy.def1.txt);
echo ${numberOfLines} # ignore the first header line

for lineNumber in {2..1001}
do
	minStress=0;
	maxStress=0;
	for fid in `ls Al99*def1.txt`
	do
		sentence=$(sed -n ${lineNumber}p Al99.NegFRho.ZeroRhoR.ZeroPhiR.eam.alloy.def1.txt);
		readStress=$(echo $sentence | awk '{print $2}');


		if [[ (1 -eq $(bc <<< "$minStress <= 0")) && (1 -eq $(bc <<< "$maxStress <= 0"))  ]]; then
			minStress=$readStress; maxStress=$readStress;
		fi
		# 						initialize 								#

		sentence=$(sed -n ${lineNumber}p ${fid} );
		readStress=$(echo $sentence | awk '{print $2}');

		# echo "${fid} -> ${lineNumber} line: $sentence. readStress = $readStress. minStress = $minStress. maxStress = ${maxStress}."
		echo "${fid} -> ${lineNumber}: readStress = $readStress. minStress = $minStress. maxStress = ${maxStress}."

		if [ 1 -eq "$(echo "${readStress} > ${maxStress}" | bc)" ]; then
			maxStress=${readStress};
		fi
		if [ 1 -eq "$(echo "${readStress} < ${minStress}" | bc)" ]; then
			minStress=${readStress};
		fi
	done
	echo "${lineNumber}: minStress = ${minStress}. maxStress = ${maxStress}.";
	echo $minStress >> $minStressFile;
	echo $maxStress >> $maxStressFile;
	minStress=0;maxStress=0;
	echo; echo;
done	
