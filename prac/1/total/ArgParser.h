#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>

class ArgParser {
	bool debug;
    bool verbose;
    int Nx;
    int Ny;
    int K1;
    int K2;
    int T;
    double tol;
    std::string input_file;
    std::string output_file;
public:
    ArgParser(): debug(false), verbose(false), Nx(0), Ny(0), K1(0), K2(0), \
				T(1), tol(0), input_file(""), output_file("") {};
    int GetNx() { return Nx; };
    int GetNy() { return Ny; };
    int GetK1() { return K1; };
    int GetK2() { return K2; };
    int GetT()  { return T;  };
    double Gettol() { return tol; };
    bool GetVerbose() { return verbose; }
	bool GetDebug() { return debug; };
    std::string GetOutput() { return output_file; };
    bool Parse(int, char**);
    void PrintHelp();
};
