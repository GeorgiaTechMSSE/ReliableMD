include ../../lss.inc

sptest: sptest.o
	$(CPP) $(ROOTDIR)/src/testgen/testgen.o sptest.o $(CXXLINKFLAGS) -o $(ROOTDIR)/bin/sptest

sptest.o: sptest.cpp
	$(CPP) $(CXXFLAGS) -c -o sptest.o sptest.cpp $(CXXOPT)

clean:
	rm -f *.o
