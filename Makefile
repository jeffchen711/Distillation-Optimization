all: distOp

distOp: distOp.o
	g++ -g distOp.o -o distOp

distOp.o: distOp.cpp distOp.h
	g++ -c -g distOp.cpp

clean:
	rm -rf *.o