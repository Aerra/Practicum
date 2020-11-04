#include "Matrices.h"
#include "kernels.h"

// Conjugate Gradient algorithm realization
void Matrices::Conjugate_Gradient(double tol, bool debug) {
    bool convergence = false;
    int counter = 1;

    double *r = new double[N1];
    for (int i = 0; i < N1; i++)
        r[i] = b[i];

    double *z = new double[N1];
    double *p = new double[N1];
    double po_0 = 0;
    double po = 0;
    double *q = new double[N1];

    do {
        kernel_VpV(rev_diag, r, z, N1);
        po = kernel_dot(N1, r, z);
        if (counter == 1) {
        //  p = z;
            kernel_axpby(N1, p, z, 0, 1);
        } else {
            //double beta = po / po_0;
            //kernel_axpby(N1, p, z, beta, 1);
            kernel_axpby(N1, p, z, po / po_0, 1);
        }
        kernel_SpMV(A, IA, JA, p, q, N1);
        double alpha = po / kernel_dot(N1, p, q);
        kernel_axpby(N1, x, p, 1, alpha);
        kernel_axpby(N1, r, q, 1, -alpha);
        po_0 = po;
        // if (po < tol || counter >= MAXITER)
		if (po < tol)
            convergence = true;
        else
            counter += 1;

		if (debug) {
    		double res_i = kernel_L2_ResidentialRate(A, IA, JA, x, b, N1);
        	std::cout << "After iter " << counter-1 << " res: " << res_i << "\n";
		}

    }
    while (!convergence);

    res = kernel_L2_ResidentialRate(A, IA, JA, x, b, N1);

	if (debug) {
    	std::cout << "X: ";
    	for (int i = 0; i < N1; i++)
    	    std::cout << " " << x[i];
    	std::cout << "\nres: " << res << "\n";
	}

    delete [] r;
    delete [] z;
    delete [] p;
    delete [] q;

    return;
}

// By already given IA and JA generate (fill) matrix A (CSR format)
void Matrices::fill(bool debug) {
    double tbeg = omp_get_wtime();
    #pragma omp parallel for    
    for (int i = 0; i < N1; i++) {
        int diag = 0;
        double sum = 0;
        for (int _j = IA[i]; _j < IA[i+1]; _j++) {
            int j = JA[_j];
            if (i == j) {
                diag = _j;
            } else {
                A[_j] = cos(i*j+i+j);
				//std::cout << A[_j] << " has ABS: " << abs(A[_j]) << " fabs " << fabs(A[_j]) << "\n";
                sum += fabs(A[_j]); // WTF is going on?? Why abs != fabs ???
            }
        }
        A[diag] = 1.234*sum;
        rev_diag[i] = 1 / A[diag]; // M^(-1)
        b[i] = sin(i);
    }
    double twcl = 0;
    twcl += omp_get_wtime() - tbeg;
	if (debug) {
    	std::cout << "Fill time: " << twcl << "\n";
	}

    return;
}

