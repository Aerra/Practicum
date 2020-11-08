CXX=g++
RM=rm -rf
CPPFLAGS=-Wall -Wextra -Wfloat-equal -O3 -fopenmp -pedantic
OBJS=ArgParser.cpp Graph.cpp kernels.cpp Matrices.cpp main.cpp

all: main

main:
	$(CXX) $(CPPFLAGS) $(OBJS) -o main 

clean:
	$(RM) main
