CXX = g++
INCL = -I/usr/include/python3.10
#INCL = -I/usr/include/python2.7
CXXFLAGS = -fPIC

all: iog.so

iog.o:iog.cpp
	g++ -c iog.cpp $(CXXFLAGS) $(INCL)

io_general.o: io_general.c	
	gcc -c io_general.c  $(CXXFLAGS) $(INCL)

io_general_class.o: io_general_class.cpp 
	g++ -c io_general_class.cpp $(CXXFLAGS) $(INCL)

iog.so:iog.o io_general.o io_general_class.o 
	g++ -o iog.so iog.o io_general_class.o io_general.o -shared $(CXXFLAGS) $(INCL)

clean:
	rm *.o *.so
