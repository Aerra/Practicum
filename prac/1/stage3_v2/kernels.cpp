#include "kernels.h"

double kernel_dot(int n, double *x, double *y) {
    double dot = 0;
    #pragma omp parallel for reduction(+:dot)
    for (int i = 0; i < n; i++)
        dot += x[i]*y[i];
    return dot;
}

void kernel_axpby(int n, double *x, double *y, double a, double b) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++)
        x[i] = a*x[i] + b*y[i];
}

void kernel_SpMV(double *A, int *IA, int *JA, double *y, double *res, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = IA[i]; j < IA[i+1]; j++) {
            int _j = JA[j];
            sum += A[j] * y[_j];
        }
        res[i] = sum;
    }
}

void kernel_VpV(double *x, double *y, double *res, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++)
        res[i] = x[i] * y[i];
}

// ||Ax-b|| = (Ax-b, Ax-b)^(1/2)
double kernel_L2_ResidentialRate(double *A, int *IA, int *JA, double *x,
									double *b, int N) {
	double *res_inter = new double[N];
	kernel_SpMV(A, IA, JA, x, res_inter, N);
	kernel_axpby(N, res_inter, b, 1, -1);
	double res = sqrt(kernel_dot(N, res_inter, res_inter));
	delete [] res_inter;
	return res;
}
