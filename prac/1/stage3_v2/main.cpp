#include <stdio.h>
#include <iostream>
#include <omp.h>
#include "kernels.h"
#include "ArgParser.h"
#include "Matrices.h"
#include "Graph.h"

void print_info(Matrices *obj, bool verbose = false, std::string output_file = "") {
	if (verbose) {
		std::cout << "IA: ";			
		for (int i = 0; i < obj->N1 + 1; i++) {
				std::cout << obj->IA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "JA: ";			
		for (int i = 0; i < obj->E1; i++) {
				std::cout << obj->JA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "N: " << obj->N1 << "\n";
		std::cout << "A: ";			
		for (int i = 0; i < obj->E1; i++) {
				std::cout << obj->A[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "b: ";			
		for (int i = 0; i < obj->N1; i++) {
				std::cout << obj->b[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "X: ";
		for (int i = 0; i < obj->N1; i++) {
				std::cout << obj->x[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "res: " << obj->res << "\n";
	}
	if (output_file != (std::string)"") {
		setlocale(LC_ALL, "rus");
		std::ofstream out_file (output_file.c_str());
		if (!out_file.is_open()) {
			std::cout << "Can't open output file " << output_file << "\n"; 
			return;
		}
		out_file << "IA: ";			
		for (int i = 0; i < obj->N1 + 1; i++) {
				out_file << obj->IA[i] << " ";
		}		
		out_file << "\n";
		out_file << "JA: ";			
		for (int i = 0; i < obj->E1; i++) {
				out_file << obj->JA[i] << " ";
		}		
		out_file << "\n";
		out_file << "N: " << obj->N1 << "\n";

		out_file << "A: ";			
		for (int i = 0; i < obj->E1; i++) {
				out_file << obj->A[i] << " ";
		}		
		out_file << "\n";
		out_file << "b: ";			
		for (int i = 0; i < obj->N1; i++) {
				out_file << obj->b[i] << " ";
		}		
		out_file << "\n";
		out_file << "X: ";
		for (int i = 0; i < obj->N1; i++) {
				out_file << obj->x[i] << " ";
		}		
		out_file << "\n";
		out_file << "res: " << obj->res << "\n";

		out_file.close();
	}
}

int main(int argc, char *argv[]) {
	// parse incoming arguments
	ArgParser Parser;
	bool continue_work = Parser.Parse(argc, argv);
	if (!continue_work) {
		return 0;
	}
	
	omp_set_num_threads(Parser.GetT());
	// GENERATE STAGE
	Graph A;
	
	A.create_graph(Parser.GetNx(), Parser.GetNy(), Parser.GetK1(),
					Parser.GetK2());
	Matrices B(A.get_N(), A.get_E());
	A.get_matrices(&B);

	// FILL STAGE
	B.fill();

	// SOLVE
	B.Conjugate_Gradient(Parser.Gettol());

	print_info(&B, Parser.GetVerbose(), Parser.GetOutput());
	return 0;
}
