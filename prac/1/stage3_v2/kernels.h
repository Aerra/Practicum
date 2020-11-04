#include <cmath>

// sum(x[i]*y[i])
double kernel_dot(int n, double *x, double *y);
// x = a*x + b*y
void kernel_axpby(int n, double *x, double *y, double a, double b);
// res = A*x
void kernel_SpMV(double *A, int *IA, int *JA, double *y, double *res, int N);
// res = x*y
void kernel_VpV(double *x, double *y, double *res, int N);
// ||Ax-b||
double kernel_L2_ResidentialRate(double *A, int *IA, int *JA, double *x, double *b, int N);
