#------------------------------------------------------------------------
# Verified solvers for parametric linear systems
#------------------------------------------------------------------------

#======== Modify the following values according to your needs ==================

#The C-XSC root directory
# CXSCDIR=/home/user/cxsc
CXSCDIR=/home/anhvt89/Documents/c++lib/cxsc-2-5-4

#The C++ compiler to use
CXX=g++

#Optimization options to use
OPTIMIZATION=-O3

#BLAS and LAPACK usage
#To use BLAS and LAPACK set this to 
#BLASLAPACK=-DCXSC_USE_BLAS -DCXSC_USE_LAPACK
#(and add Include-Paths or Compiler Options for your BLAS/LAPACK
#library if necessary)
#otherwise leave blank
BLASLAPACK=-DCXSC_USE_BLAS -DCXSC_USE_LAPACK

#Activate OpenMP by setting this flag to the appropriate
#OpenMP-option for your compiler (GCC: -fopenmp, INTEL: -openmp)
OPENMP=-fopenmp

#Command for generating static libraries
LIBTOOL=ar rcs


#======== The following options are only necessary for compiling the examples =====

#Link options to link against C-XSC
#The -Wl,-R option sets the search path for a dynamic
#C-XSC library. If you use a static C-XSC library this
#option is not necessary. Depending on your system/compiler
#a different compiler option might be necessary.
CXSCLINK=-Wl,-R$(CXSCDIR)/lib -L$(CXSCDIR)/lib -lcxsc

#Further options to use for linking
#If you have activated BLAS and LAPACK usage
#you must link to appropriate BLAS and LAPACK libraries
LINKOPTS=-llapack -lblas 



#======== Do not modify beyond this line unless you know what you do ===========

CXXOPTS=$(OPTIMIZATION) $(BLASLAPACK) $(OPENMP)

SRCDIR=./src
LIBDIR=./lib
BINDIR=./bin
INCDIR=./include

all: lib examples

$(SRCDIR)/parlinsys.o: $(SRCDIR)/parlinsys.cpp $(INCDIR)/parlinsys.hpp
	$(CXX) $(CXXOPTS) -c -I. -I$(CXSCDIR)/include -I$(INCDIR) $(SRCDIR)/parlinsys.cpp -o $(SRCDIR)/parlinsys.o

$(SRCDIR)/symlinsys.o: $(SRCDIR)/symlinsys.cpp $(INCDIR)/symlinsys.hpp $(INCDIR)/parlinsys.hpp
	$(CXX) $(CXXOPTS) -c -I. -I$(CXSCDIR)/include -I$(INCDIR) $(SRCDIR)/symlinsys.cpp -o $(SRCDIR)/symlinsys.o

$(SRCDIR)/parlinsysnp.o: $(SRCDIR)/parlinsysnp.cpp $(INCDIR)/parlinsys.hpp
	$(CXX) $(CXXOPTS) -c -I. -I$(CXSCDIR)/include -I$(INCDIR) $(SRCDIR)/parlinsysnp.cpp -o $(SRCDIR)/parlinsysnp.o


lib: $(SRCDIR)/parlinsys.o $(SRCDIR)/parlinsysnp.o $(SRCDIR)/symlinsys.o 
	$(LIBTOOL) $(LIBDIR)/libparlinsys.a $(SRCDIR)/parlinsys.o $(SRCDIR)/parlinsysnp.o $(SRCDIR)/symlinsys.o 

examples: $(BINDIR)/parlinsys_example $(BINDIR)/parlinsys_example2 $(BINDIR)/symlinsys_example

$(BINDIR)/parlinsys_example: $(SRCDIR)/example.cpp lib
	$(CXX) $(CXXOPTS) $(SRCDIR)/example.cpp -o $(BINDIR)/parlinsys_example -I$(INCDIR) -I$(CXSCDIR)/include -L$(LIBDIR) $(LINKOPTS) $(CXSCLINK) -lparlinsys

$(BINDIR)/parlinsys_example2: $(SRCDIR)/example2.cpp lib
	$(CXX) $(CXXOPTS) $(SRCDIR)/example2.cpp -o $(BINDIR)/parlinsys_example2 -I$(INCDIR) -I$(CXSCDIR)/include -L$(LIBDIR) $(LINKOPTS) $(CXSCLINK) -lparlinsys

$(BINDIR)/symlinsys_example: $(SRCDIR)/example_sym.cpp lib
	$(CXX) $(CXXOPTS) $(SRCDIR)/example_sym.cpp -o $(BINDIR)/symlinsys_example -I$(INCDIR) -I$(CXSCDIR)/include -L$(LIBDIR) $(LINKOPTS) $(CXSCLINK) -lparlinsys

clean:
	rm -f $(BINDIR)/parlinsys_example $(BINDIR)/parlinsys_example2 $(BINDIR)/symlinsys_example
	rm -f $(LIBDIR)/libparlinsys.a
	rm -f $(SRCDIR)/parlinsys.o $(SRCDIR)/parlinsysnp.o $(SRCDIR)/symlinsys.o $(SRCDIR)/*~

