#------------------------------------------------------------------------
# Verified solvers for parametric linear systems
#------------------------------------------------------------------------

#======== Modify the following values according to your needs ==================

#The C-XSC root directory
# CXSCDIR=/home/user/cxsc
CXSCDIR=/home/anhvt89/Documents/c++lib/cxsc-2-5-4

#The SuiteSparse root directory
# SUITESPARSEDIR=/home/user/SuiteSparse
SUITESPARSEDIR=/home/anhvt89/Documents/c++lib/cxsc-2-5-4/examples/SparseLinSys

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
UMFPACKDIR = $(SUITESPARSEDIR)/UMFPACK
AMDDIR = $(SUITESPARSEDIR)/AMD
CAMDDIR = $(SUITESPARSEDIR)/CAMD
UFCONFIGDIR = $(SUITESPARSEDIR)/UFConfig
CHOLMODDIR = $(SUITESPARSEDIR)/CHOLMOD
COLAMDDIR=$(SUITESPARSEDIR)/COLAMD
CCOLAMDDIR=$(SUITESPARSEDIR)/CCOLAMD
METISDIR=$(SUITESPARSEDIR)/metis-4.0
SPCONFIGDIR=$(SUITESPARSEDIR)/SuiteSparse_config

LIBS=-L$(UMFPACKDIR)/Lib -L$(AMDDIR)/Lib -L$(CHOLMODDIR)/Lib -L$(COLAMDDIR)/Lib  -L$(CCOLAMDDIR)/Lib -L$(CAMDDIR)/Lib -L$(METISDIR) -L$(SPCONFIGDIR) -lumfpack -lcholmod -lamd -lcamd -lcolamd -lccolamd -lmetis -lsuitesparseconfig -lcxsc -lm -lrt

CXXOPTS=$(OPTIMIZATION) $(BLASLAPACK) $(OPENMP)

SRCDIR=./src
LIBDIR=./lib
BINDIR=./bin
INCDIR=./include

all: lib examples

$(SRCDIR)/sparselinsyskraw.o: $(SRCDIR)/sparselinsyskraw.cpp $(INCDIR)/sparselinsys.hpp $(SRCDIR)/cxsc_umfpack.hpp $(SRCDIR)/trisolve.hpp $(SRCDIR)/utility.hpp $(SRCDIR)/cxsc_cholmod.hpp
	$(CXX) $(CXXOPTS) -c -I$(SRCDIR) -I$(INCDIR) -I$(CXSCDIR)/include -I$(AMDDIR)/Include -I$(UMFPACKDIR)/Include -I$(CHOLMODDIR)/Include -I$(UFCONFIGDIR) -I$(SPCONFIGDIR) $(SRCDIR)/sparselinsyskraw.cpp -o $(SRCDIR)/sparselinsyskraw.o

$(SRCDIR)/trisolve.o: $(SRCDIR)/trisolve.hpp $(SRCDIR)/trisolve.cpp $(SRCDIR)/QRdecomp.hpp $(SRCDIR)/utility.hpp
	$(CXX) $(CXXOPTS) -c -I$(SRCDIR) -I$(INCDIR) -I$(CXSCDIR)/include -I$(AMDDIR)/Include -I$(UMFPACKDIR)/Include -I$(CHOLMODDIR)/Include -I$(UFCONFIGDIR) -I$(SPCONFIGDIR) $(SRCDIR)/trisolve.cpp -o $(SRCDIR)/trisolve.o

$(SRCDIR)/sparselinsyssing.o: $(SRCDIR)/sparselinsyssing.cpp $(INCDIR)/sparselinsys.hpp $(SRCDIR)/cxsc_umfpack.hpp $(SRCDIR)/trisolve.hpp $(SRCDIR)/utility.hpp $(SRCDIR)/cxsc_cholmod.hpp
	$(CXX) $(CXXOPTS) -c -I$(SRCDIR) -I$(INCDIR) -I$(CXSCDIR)/include -I$(AMDDIR)/Include -I$(UMFPACKDIR)/Include -I$(CHOLMODDIR)/Include -I$(UFCONFIGDIR) -I$(SPCONFIGDIR) $(SRCDIR)/sparselinsyssing.cpp -o $(SRCDIR)/sparselinsyssing.o -DCXSC_FAST_COMPLEX_OPERATIONS

$(SRCDIR)/utility.o: $(SRCDIR)/utility.cpp $(SRCDIR)/utility.hpp
	$(CXX) $(CXXOPTS) -c -I$(SRCDIR) -I$(INCDIR) -I$(CXSCDIR)/include -I$(AMDDIR)/Include -I$(CHOLMODDIR)/Include -I$(UFCONFIGDIR) -I$(SPCONFIGDIR) $(SRCDIR)/utility.cpp -o $(SRCDIR)/utility.o -DCXSC_FAST_COMPLEX_OPERATIONS

$(SRCDIR)/cxsc_cholmod.o: $(SRCDIR)/cxsc_cholmod.cpp $(SRCDIR)/cxsc_cholmod.hpp
	$(CXX) $(CXXOPTS) -c -I$(SRCDIR) -I$(INCDIR) -I$(CXSCDIR)/include -I$(AMDDIR)/Include -I$(CHOLMODDIR)/Include -I$(UFCONFIGDIR) -I$(SPCONFIGDIR) $(SRCDIR)/cxsc_cholmod.cpp -o $(SRCDIR)/cxsc_cholmod.o -DCXSC_FAST_COMPLEX_OPERATIONS

$(SRCDIR)/cxsc_umfpack.o: $(SRCDIR)/cxsc_umfpack.cpp $(SRCDIR)/cxsc_umfpack.hpp
	$(CXX) $(CXXOPTS) -c -I$(SRCDIR) -I$(INCDIR) -I$(CXSCDIR)/include -I$(AMDDIR)/Include -I$(UMFPACKDIR)/Include -I$(CHOLMODDIR)/Include -I$(UFCONFIGDIR) -I$(SPCONFIGDIR) $(SRCDIR)/cxsc_umfpack.cpp -o $(SRCDIR)/cxsc_umfpack.o -DCXSC_FAST_COMPLEX_OPERATIONS

lib: $(SRCDIR)/sparselinsyskraw.o $(SRCDIR)/sparselinsyssing.o $(SRCDIR)/trisolve.o $(SRCDIR)/utility.o $(SRCDIR)/cxsc_cholmod.o $(SRCDIR)/cxsc_umfpack.o
	$(LIBTOOL) $(LIBDIR)/libsparselinsys.a $(SRCDIR)/sparselinsyskraw.o $(SRCDIR)/sparselinsyssing.o $(SRCDIR)/trisolve.o $(SRCDIR)/utility.o $(SRCDIR)/cxsc_cholmod.o $(SRCDIR)/cxsc_umfpack.o 

examples: $(BINDIR)/sparselinsys_simple_example $(BINDIR)/sparselinsys_mm_example

$(BINDIR)/sparselinsys_simple_example: $(SRCDIR)/simple_example.cpp $(INCDIR)/sparselinsys.hpp lib
	$(CXX) $(CXXOPTS) $(SRCDIR)/simple_example.cpp -o $(BINDIR)/sparselinsys_simple_example -I$(INCDIR) -I$(CXSCDIR)/include $(LINKOPTS) -L$(LIBDIR) $(CXSCLINK) -lsparselinsys $(LIBS)

$(BINDIR)/sparselinsys_mm_example: $(SRCDIR)/mm_example.cpp $(INCDIR)/sparselinsys.hpp lib
	$(CXX) $(CXXOPTS) $(SRCDIR)/mm_example.cpp -o $(BINDIR)/sparselinsys_mm_example -I$(INCDIR) -I$(CXSCDIR)/include $(LINKOPTS) -L$(LIBDIR) $(CXSCLINK) -lsparselinsys $(LIBS)

clean:
	rm -f $(BINDIR)/sparselinsys_mm_example $(BINDIR)/sparselinsys_simple_example
	rm -f $(LIBDIR)/libsparselinsys.a
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*~

