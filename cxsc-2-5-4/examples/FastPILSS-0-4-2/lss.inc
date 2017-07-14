#-------------------------------------------------------------------------------
# C-XSC install directory
#-------------------------------------------------------------------------------
# PREFIX=/home/user/cxsc
PREFIX=/home/anhvt89/Documents/c++lib/cxsc-2-5-4

#-------------------------------------------------------------------------------
# The rootdir of the programms (=the directory containing /src) 
#-------------------------------------------------------------------------------
# ROOTDIR=/home/user/FastPILSS
ROOTDIR=/home/anhvt89/Documents/c++lib/cxsc-2-5-4/FastPILSS-0-4-2/

#-------------------------------------------------------------------------------
# C++ compiler for serial programms an MPI-C++ compiler for parallel programms
#-------------------------------------------------------------------------------
CPP=g++ 
MPICPP=mpic++

#-------------------------------------------------------------------------------
# Options for optimization and inlining (should always be activated)
#-------------------------------------------------------------------------------
OPTIMIZATION=-O3 

#-------------------------------------------------------------------------------
# Flag for using OpenMP
#-------------------------------------------------------------------------------
OPENMP=-fopenmp

#-------------------------------------------------------------------------------
# BLAS and LAPACK link commands
# (order of the linked libraries might be important!)
#-------------------------------------------------------------------------------
BLASVERSION=-llapack -lblas

#-------------------------------------------------------------------------------
# Link flags for ScaLAPACK and BLACS
#-------------------------------------------------------------------------------
SCALAPACKLINK=-lscalapack -lblacs

#-------------------------------------------------------------------------------
# Link flags for fortran (if necessary)
#-------------------------------------------------------------------------------
FORTRAN=-lgfortran

#-------------------------------------------------------------------------------
# Naming convention when using Fortran functions from C/C++
# Possible choices:
# Unchanged: Function names are the same in C/C++ and Fortran
# Add_: A _ is added to the function name ( for example dgemm_ )
# Uppercase: Function names are upper case ( for example DGEMM )
# Change this setting if you get linker errors concerning ScaLAPACK/BLACS 
# functions
#-------------------------------------------------------------------------------
FORTRANNAMING=Add_

#-------------------------------------------------------------------------------
# Comment this line or leave it blank if you compiled C-XSC into a static 
# library
#-------------------------------------------------------------------------------
CXXRPATH=-Wl,-R$(PREFIX)/lib                          

#-------------------------------------------------------------------------------
# Command for creating a static library on your system
# For Intel compiler using -ipo option use xiar instead of ar
#-------------------------------------------------------------------------------
LIBTOOL=ar rcs

#-------------------------------------------------------------------------------
# The make command you want to use (gmake or make in most cases)
#-------------------------------------------------------------------------------
MAKETOOL=make



#-------------------------------------------------------------------------------
# Do not edit below this line unless you know what you do!
#-------------------------------------------------------------------------------

CXXOPTS=$(OPTIMIZATION) $(OPENMP) -D$(FORTRANNAMING)
CXXINC=-I$(PREFIX)/include  -I$(ROOTDIR)/src/pdotk -I$(ROOTDIR)/src/testgen -I$(ROOTDIR)/src/lss
CXSCLIBOPTS=-L$(PREFIX)/lib $(CXXRPATH)  -lcxsc 

CXXFLAGS=$(CXXINC) $(CXXOPTS)              
CXXLINKFLAGS=$(CXSCLIBOPTS) $(BLASVERSION) $(FORTRAN) $(OPENMP)
CXXMPILINKFLAGS=$(CXSCLIBOPTS) -lm -lpthread $(SCALAPACKLINK) $(BLASVERSION) $(OPENMP) $(FORTRAN)

