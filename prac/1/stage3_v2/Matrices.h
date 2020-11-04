//This header describe different types of matrices that we use for CSR format
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
    };
    ~Matrices() {
        delete[] IA;
        delete[] JA;
        delete[] A;
        delete[] b;
        delete[] x;
        delete[] rev_diag;
    };
    void fill(bool=false); // debug = false
    void Conjugate_Gradient(double, bool=false); // tol, debug = false
};

