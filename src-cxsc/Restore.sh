#!/bin/bash
# restore or store a file


file=pair_eam_rad.h
# version=.bak.11Feb16.B
version=.bak.1Feb16

cp ${file}${version} ${file} # restore from back up
echo ""
echo "${file}${version}"
echo "-> ${file} "
echo ""

# cp ${file} ${file}${version} # save to restore
# echo ""
# echo "${file} "
# echo "->${file}${version}"
# echo ""