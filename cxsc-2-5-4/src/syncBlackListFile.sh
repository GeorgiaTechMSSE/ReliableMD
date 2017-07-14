#!/bin/bash

fileList=$(cat blacklist.txt);
numberOfLines=$(wc -l < blacklist.txt);
srcDirectory=/home/anhvt89/Documents/tar/cxsc-2-5-4/src/
echo "${fileList}"
echo "number of line: ${numberOfLines}"
echo ""
for lineNumber in $(seq $numberOfLines)
do
	fileName=$(sed -n ${lineNumber}p blacklist.txt);
	rm ${fileName} || true;
	cp ${srcDirectory}/${fileName} . 
	echo "lineNumber ${lineNumber}: ${fileName} - copy = done!";
done
