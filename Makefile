OPTIONS = -O2 -Wall -Wextra -pedantic -ansi -DNDEBUG
LIBS = -lgsl -lgslcblas -ldnest3 -lboost_thread-mt -lboost_system

default:
	g++ $(OPTIONS) -c *.cpp
	g++ -o main *.o $(LIBS)
	rm -f main.o
	ar rcs libstarfield.a *.o
	rm -f *.o

