include ../../lss.inc

plss: plss_ex.o
	$(MPICPP) $(ROOTDIR)/src/examples/plss_ex.o $(ROOTDIR)/src/testgen/testgen.o $(ROOTDIR)/lib/libcxscplss.a $(CXXMPILINKFLAGS) -o $(ROOTDIR)/bin/plss 

pilss: pilss_ex.o
	$(MPICPP) $(ROOTDIR)/src/examples/pilss_ex.o $(ROOTDIR)/src/testgen/testgen.o $(ROOTDIR)/lib/libcxscplss.a $(CXXMPILINKFLAGS) -o $(ROOTDIR)/bin/pilss 

pclss: pclss_ex.o
	$(MPICPP) $(ROOTDIR)/src/examples/pclss_ex.o $(ROOTDIR)/src/testgen/testgen.o $(ROOTDIR)/lib/libcxscplss.a $(CXXMPILINKFLAGS) -o $(ROOTDIR)/bin/pclss 

pcilss: pcilss_ex.o
	$(MPICPP) $(ROOTDIR)/src/examples/pcilss_ex.o $(ROOTDIR)/src/testgen/testgen.o $(ROOTDIR)/lib/libcxscplss.a $(CXXMPILINKFLAGS) -o $(ROOTDIR)/bin/pcilss 

plss_ex.o: plss.cpp $(ROOTDIR)/src/testgen/testgen.hpp
	$(MPICPP) $(CXXFLAGS) -I$(ROOTDIR)/include -c -o plss_ex.o plss.cpp $(CXXOPT) 
pclss_ex.o: pclss.cpp $(ROOTDIR)/src/testgen/testgen.hpp
	$(MPICPP) $(CXXFLAGS) -I$(ROOTDIR)/include -c -o pclss_ex.o pclss.cpp $(CXXOPT)
pilss_ex.o: pilss.cpp $(ROOTDIR)/src/testgen/testgen.hpp
	$(MPICPP) $(CXXFLAGS) -I$(ROOTDIR)/include -c -o pilss_ex.o pilss.cpp $(CXXOPT)
pcilss_ex.o: pcilss.cpp $(ROOTDIR)/src/testgen/testgen.hpp
	$(MPICPP) $(CXXFLAGS) -I$(ROOTDIR)/include -c -o pcilss_ex.o pcilss.cpp $(CXXOPT)
	
clean:
	rm -f *.o
  	
