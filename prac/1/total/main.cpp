#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include "kernels.h"
#include "ArgParser.h"
#include "Matrices.h"
#include "Graph.h"

void print_info(Matrices *obj, bool verbose = false, \
			   	std::string output_file = "") {
	if (verbose) {
		std::cout << "\nIA: ";			
		for (int i = 0; i < obj->N1 + 1; i++) {
				std::cout << obj->IA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "JA: ";			
		for (int i = 0; i < obj->E1; i++) {
				std::cout << obj->JA[i] << " ";
		}		
		std::cout << "\nE: " << obj->E1;
		std::cout << "\nN: " << obj->N1;
		std::cout << "\nA: ";			
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

void print_generate(int N, double time) {
	std::cout << "\n+-------------------------------------+\n";
	std::cout << "|                 GENERATE            |\n";
	std::cout << "+-------------------------------------+\n";
	std::cout << "| N (count vertices) |" << std::setw(16) << N << "|\n";
	std::cout << "+-------------------------------------+\n";
	std::cout << "| Wasted Time        |" << std::setw(15) << time << "s|\n";
	std::cout << "+-------------------------------------+\n";
}

void print_fill(double time) {
	std::cout << "\n+-------------------------------------+\n";
	std::cout << "|               FILL                  |\n";
	std::cout << "+-------------------------------------+\n";
	std::cout << "| Wasted Time        |" << std::setw(15) << time << "s|\n";
	std::cout << "+-------------------------------------+\n";
}

void print_solve(int iter, double res, double time, double dot_time, \
			   	double axpby_time, double SpMV_time, double VpV_time) {
	std::cout << "\n+-----------------------------------------+\n";
	std::cout << "|                     SOLVE               |\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| Iterations             |" << std::setw(16) << iter << "|\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| ||Ax-b||               |" << std::setw(16) << res << "|\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| Wasted Time            |" << std::setw(15) << time << "s|\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| where:                                  |\n";
	std::cout << "| 2*Iterations for dot   |" << std::setw(15) << dot_time 
				<< "s|\n";
	std::cout << "+-----------------------------------------+\n";
	double av_dot = dot_time / (2 * iter);
	std::cout << "| Average dot time       |" << std::setw(15)
				<< av_dot << "s|\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| 3*Iterations for axpby |" << std::setw(15) << axpby_time
				<< "s|\n";
	std::cout << "+-----------------------------------------+\n";
	double av_axpby = axpby_time / (3 * iter);
	std::cout << "| Average axpby time     |" << std::setw(15)
				<< av_axpby << "s|\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| 1*Iterations for SpMV  |" << std::setw(15) << SpMV_time
				<< "s|\n";
	std::cout << "+-----------------------------------------+\n";
	double av_SpMV = SpMV_time / iter;
	std::cout << "| Average SpMV time      |" << std::setw(15)
				<< av_SpMV << "s|\n";
	std::cout << "+-----------------------------------------+\n";
	std::cout << "| 1*Iterations for VpV   |" << std::setw(15) << VpV_time
				<< "s|\n";
	std::cout << "+-----------------------------------------+\n";
	double av_VpV = VpV_time / iter;
	std::cout << "| Average VpV time       |" << std::setw(15)
				<< av_VpV << "s|\n";
	std::cout << "+-----------------------------------------+\n";
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
	// start timing
	double generate_time = 0;
	Graph A((Parser.GetNx() + 1)*(Parser.GetNy() + 1));
	
	A.create_graph(Parser.GetNx(), Parser.GetNy(), Parser.GetK1(),
					Parser.GetK2(), generate_time, Parser.GetDebug());
	Matrices B(A.get_N(), A.get_E());
	A.get_matrices(&B, generate_time, Parser.GetT(), Parser.GetDebug());
	print_generate(A.get_N(), generate_time);

	// FILL STAGE
	double fill_time = 0;
	B.fill(fill_time, Parser.GetDebug());
	print_fill(fill_time);

	// SOLVE
	double dot_time = 0;
	double axpby_time = 0;
	double SpMV_time = 0;
	double VpV_time = 0;

	double solve_time = 0;
	B.Conjugate_Gradient(Parser.Gettol(), solve_time, dot_time, axpby_time, 
					SpMV_time, VpV_time, Parser.GetDebug());
	print_solve(B.iterations, B.res, solve_time, dot_time, axpby_time, 
				SpMV_time, VpV_time);

	print_info(&B, Parser.GetVerbose(), Parser.GetOutput());
	return 0;
}
