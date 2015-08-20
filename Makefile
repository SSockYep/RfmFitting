all: Source3

Source3: Source3.cpp RFMCalc.o
	g++ -o Source3 Source3.cpp RFMCalc.o

RFMCalc.o: RFMCalc.cpp
	g++ -c RFMCalc.cpp

clean:
	rm *.o Source3


