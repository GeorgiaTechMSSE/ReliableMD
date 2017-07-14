include ../../lss.inc

lib: matmul_par.o cxsc_pblas.o matinv_aprx_par.o mpi_comm.o templplss.o
	$(LIBTOOL) $(ROOTDIR)/lib/libcxscplss.a matmul_par.o cxsc_pblas.o matinv_aprx_par.o mpi_comm.o templplss.o 

matinv_aprx_par.o: matinv_aprx_par.cpp matinv_aprx_par.hpp
	$(MPICPP)  $(CXXFLAGS) -c -o matinv_aprx_par.o matinv_aprx_par.cpp $(CXXOPT)
matmul_par.o: matmul_par.cpp matmul_par.hpp cxsc_mpicomm.hpp
	$(MPICPP)  $(CXXFLAGS) -c -o matmul_par.o matmul_par.cpp $(CXXOPT)
mpi_comm.o: cxsc_mpicomm.cpp cxsc_mpicomm.hpp
	$(MPICPP)  $(CXXFLAGS) -c -o mpi_comm.o cxsc_mpicomm.cpp $(CXXOPT)         
templplss.o: templplss.cpp templplss.hpp matinv_aprx_par.hpp matmul_par.hpp cxsc_pblas.hpp matmul_par.hpp cxsc_mpicomm.hpp
	$(MPICPP)  $(CXXFLAGS) -c -o templplss.o templplss.cpp $(CXXOPT)
cxsc_pblas.o: cxsc_pblas.cpp cxsc_pblas.hpp
	$(MPICPP)  $(CXXFLAGS) -c -o cxsc_pblas.o cxsc_pblas.cpp $(CXXOPT)

clean:
	rm -f *.o









