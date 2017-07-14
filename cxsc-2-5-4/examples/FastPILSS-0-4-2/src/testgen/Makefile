include ../../lss.inc

tester: tester.o
	$(CPP) tester.o testgen.o $(CXXLINKFLAGS) -o $(ROOTDIR)/bin/tester

testgen.o: testgen.cpp testgen.hpp
	$(CPP)  $(CXXFLAGS) -c -o testgen.o testgen.cpp $(CXXOPT)

tester.o: test.cpp
	$(CPP)  $(CXXFLAGS) -c -o tester.o test.cpp $(CXXOPT)
	
clean:
	rm -f *.o
