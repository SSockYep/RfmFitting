all: Source1 Source2

Source1: Source1.cpp RFMCalc.o
	g++ -o Source1 Source1.cpp RFMCalc.o

Source2: Source2.cpp RFMCalc.o
	g++ -o Source2 Source2.cpp RFMCalc.o

RFMCalc.o: RFMCalc.cpp
	g++ -c RFMCalc.cpp

clean:
	rm *.o Source1 Source2


