#!/bin/bash

cd /home/anhvt89/Documents/lammps-28Jun14/src-midpoint-radius
make openmpi

cd /home/anhvt89/Documents/lammps-28Jun14/src-cxsc
make openmpi

cd /home/anhvt89/Documents/lammps-28Jun14/src-original
make openmpi

cd /home/anhvt89/Documents/lammps-28Jun14/src-total-uncertainty-abs
make openmpi

cd /home/anhvt89/Documents/lammps-28Jun14/src-total-uncertainty-noabs
make openmpi

cd /home/anhvt89/Documents/lammps-28Jun14/src-ublb
make openmpi