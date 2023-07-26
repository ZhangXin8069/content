CXX = g++
INCL = -I/dssg/opt/icelake/linux-centos8-icelake/gcc-8.3.1/python-3.9.12-omx7ttx5qbayml44yyhwifxuooiebvh7/include
#INCL = -I/usr/include/python3.6m
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
