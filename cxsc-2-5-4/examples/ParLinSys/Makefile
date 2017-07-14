#------------------------------------------------------------------------
# Makefile for the module parlinsys and corresponding sample programs of 
# the parametric Supplement to "C++ Toolbox for Verified  Computing I".
#------------------------------------------------------------------------

#======== modify the following values ==================

#the C-XSC root directory
CXSCDIR  = /home/user/cxsc

# which C++ compiler? 
# Possible values: gnu         for the GNU C/C++ Compiler (gcc)
#                  sun_forte   for the SUN Forte Compiler (OS SUN Solaris)
#                  intel       for the INTEL C/C++ Compiler (OS Linux)
COMPILER=gnu

#========= the following commands should work on most Unix systems ========

# settings for SUN Forte C++ Compiler --------------------------------
ifeq ($(COMPILER),sun_forte)
CXX=CC#                          name of the C++ compiler
CXXOPTS=#                        optional flags forwarded to the C++ compiler
LIBS=-lcxsc -lm9x#     linking with libraries

else
# settings for Intel C++ Compiler ------------------------------------
ifeq ($(COMPILER),intel)
CXX=icc#                         name of the C++ compiler
CXXOPTS=#                        optional flags forwarded to the C++ compiler
LIBS=-lcxsc -lm#       linking with libraries

else
# settings for GNU - Compiler ----------------------------------------
CXX=g++#                         name of the C++ compiler
CXXOPTS=-Wall -Winline#          optional flags forwarded to the C++ compiler
LIBS=-Wl,-R$(CXSCDIR)/lib -lcxsc -lm#       linking with libraries

endif
endif

# Settings for all systems -------------------------------------------

# Suffix rules
#-------------
.SUFFIXES:
.SUFFIXES: .cpp .o .a 

all: parlinsys parlinsys_exf

parlinsys: parlinsys.cpp
	$(CXX) $(CXXOPTS) -c -I. -I$(CXSCDIR)/include parlinsys.cpp

parlinsys_exf: parlinsys.o parlinsys_exf.cpp
	$(CXX) $(CXXOPTS) -o parlinsys_exf -I. -I$(CXSCDIR)/include -L$(CXSCDIR)/lib parlinsys_exf.cpp parlinsys.o $(LIBS)

clean:
	rm -f parlinsys_exf
	rm -f parlinsys.o

