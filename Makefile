CC=g++
LDLIBS= -L./. -lRFMCalc_Release -L/home/gsl -lgsl -lgslcblas

OBJECTS=Source1.o libRFMCalc_Release.a

Source1:    $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDLIBS) 
	chmod 755 $@

Source1.o :	Source1.cpp RFMCalcDll.h
	$(CC) -c -O3 Source1.cpp -lm -ltiff
	
libRFMCalc_Release.a: Mat.o RFMCalcDll.o
	ar rcs libRFMCalc_Release.a Mat.o RFMCalcDll.o
	
Mat.o: Mat.cpp
	g++ -c -o Mat.o Mat.cpp
RFMCalcDll.o: RFMCalcDll.cpp
	g++ -c -o RFMCalcDll.o RFMCalcDll.cpp

clean:
	rm -f Source1 
	rm -f *.o
	rm -f *.a