#include <cmath>
#include <omp.h>

// sum(x[i]*y[i])
double kernel_dot(int n, double *x, double *y, double &twcl);
// x = a*x + b*y
void kernel_axpby(int n, double *x, double *y, double a, double b,\
			   	double &twcl);
// res = A*y
void kernel_SpMV(double *A, int *IA, int *JA, double *y, double *res, int N, \
				double &twcl);
// res = x*y
void kernel_VpV(double *x, double *y, double *res, int N, double &twcl);
// ||Ax-b||
double kernel_L2_ResidentialRate(double *A, int *IA, int *JA, double *x, \
				double *b, int N);
