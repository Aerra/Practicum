Program solves A1 variant

Compile:
g++ -Wall -O3 -fopenmp -pedantic ArgParser.cpp Graph.cpp kernels.cpp Matrices.cpp main.cpp -o main

Help: ./main -h

Launch:
1) Arguments from command line
./main --Nx=10 --Ny=10 --K1=2 --K2=3 --T=4 --tol=0.001 --output="example/output.txt"

2) Arguments from input file
./main --input="examples/input.txt"
