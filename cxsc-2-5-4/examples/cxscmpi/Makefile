CXSCDIR=/home/anhvt89/Documents/c++lib/cxsc-2-5-4
CXSCLINK=-L$(CXSCDIR) -Wl,-R$(CXSCDIR) # dynamic linking
SRCDIR=$(CXSCDIR)/src
LIBDIR=$(CXSCDIR)/lib
BINDIR=$(CXSCDIR)/bin
INCDIR=$(CXSCDIR)/include

default: build

build: cxsc_mpi_example.o cxsc_mpicomm.o
	mpic++ \
		-W -Wall -Wextra -pedantic -std=c++0x \
		-I$(CXSCDIR) -I$(CXSCDIR)/include -I$(INCDIR) -I$(SRCDIR) \
		$(CXSCLINK)
		-lcxsc
		-lmpi \
		cxsc_mpi_example.o \
		cxsc_mpicomm.o \
		-o cxsc_mpi_example

cxsc_mpicomm.o: cxsc_mpicomm.cpp cxsc_mpicomm.hpp
	mpic++ \
		-W -Wall -Wextra -pedantic -std=c++0x \
		-I$(CXSCDIR) -I$(CXSCDIR)/include -I$(INCDIR) -I$(SRCDIR) \
		$(CXSCLINK)
		-lcxsc
		-lmpi \
		-c cxsc_mpicomm.cpp \
		-o cxsc_mpicomm.o

cxsc_mpi_example.o: cxsc_mpi_example.cpp
	mpic++ \
		-W -Wall -Wextra -pedantic -std=c++0x \
		-I$(CXSCDIR) -I$(CXSCDIR)/include -I$(INCDIR) -I$(SRCDIR) \
		$(CXSCLINK)
		-lcxsc
		-lmpi \
		-c cxsc_mpi_example.cpp \
		-o cxsc_mpi_example.o
custom: cxsc_mpicomm.cpp cxsc_mpicomm.hpp cxsc_mpi_example.cpp 
	mpic++ \
		-W -Wall -Wextra -pedantic -std=c++0x \
		-I$(CXSCDIR) -I$(CXSCDIR)/include -I$(INCDIR) -I$(SRCDIR) \
		$(CXSCLINK)
		-lcxsc
		-lmpi \
		-c cxsc_mpicomm.cpp cxsc_mpicomm.hpp cxsc_mpi_example.cpp 
		-o cxsc_mpi_example_custom
clean:
	rm -fr ./*.o

.PHONY: clean
