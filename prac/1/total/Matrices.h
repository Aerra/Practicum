//This header describe different types of matrices that we use for CSR format
// pragma once because include "Martices.h" you can find not only in one file
#pragma once
#include <iostream>
#include <omp.h>
#include <cmath>

struct Matrices {
    int *IA; // rowptr for CSR
    int *JA; // colptr for CSR
    double *A; // matrix coefficients for CSR
    double *b; // right side vector 
    double *x; // solution vector
    double *rev_diag; // M^(-1) diagonal matrix for CG 
    double res; // L2-residual rate for solution
	int iterations; // count of iterations for solver
    int N1; // matrices size (N1xN1)
    int E1; // JA length
    Matrices(int N, int E) {
        IA = new int[N+1];
        JA = new int[E];
        A = new double[E];
        b = new double[N];
        x = new double[N];
        rev_diag = new double[N];
        N1 = N;
        E1 = E;
        res = 0;
		iterations = 1;
    };
    ~Matrices() {
        delete[] IA;
        delete[] JA;
        delete[] A;
        delete[] b;
        delete[] x;
        delete[] rev_diag;
    };
    void fill(double &, bool=false); // debug = false double & for timings
	// Args: tol, twcl (timing), timing for dot, timing for axpby,
	// timing for SpMV, timing for VpV, debug = false
    void Conjugate_Gradient(double, double &, double &, double &, double &,
							double &, bool=false);
};

