include lss.inc

all: library examples tester

examples: plss pilss pclss pcilss

library: 
	( cd src ; cd plss ; $(MAKETOOL) lib )

testgen:
	( cd src ; cd testgen ; $(MAKETOOL) testgen.o )
        
sptest: testgen 
	( cd src ; cd sptest ; $(MAKETOOL) sptest )

tester: testgen  
	( cd src ; cd testgen ; $(MAKETOOL) tester )

plss: testgen library
	( cd src ; cd examples ; $(MAKETOOL) plss )

pclss: testgen library
	( cd src ; cd examples ; $(MAKETOOL) pclss )
	
pilss: testgen library
	( cd src ; cd examples ; $(MAKETOOL) pilss )

pcilss: testgen library
	( cd src ; cd examples ; $(MAKETOOL) pcilss )


clean:
	( cd bin ; rm -f sptest tester plss pilss pclss pcilss )
	( cd src ; cd examples ; $(MAKETOOL) clean )
	( cd src ; cd sptest ; $(MAKETOOL) clean )
	( cd src ; cd testgen ; $(MAKETOOL) clean )
	( cd src ; cd plss ; $(MAKETOOL) clean )
	( cd lib ; rm -f *.a *.so)











