#include "Matrices.h"
#include "kernels.h"

// Conjugate Gradient algorithm realization
void Matrices::Conjugate_Gradient(double tol, double &twcl, double &dot_twcl, \
								double &axpby_twcl, double &SpMV_twcl, \
								double &VpV_twcl, bool debug) {
    bool convergence = false;

    double *r = new double[N1];
    for (int i = 0; i < N1; i++)
        r[i] = b[i];

    double *z = new double[N1];
    double *p = new double[N1];
    double po_0 = 0;
    double po = 0;
    double *q = new double[N1];

    double tbeg = omp_get_wtime();
    do {
        kernel_VpV(rev_diag, r, z, N1, VpV_twcl);
        po = kernel_dot(N1, r, z, dot_twcl);
        if (iterations == 1) {
        //  p = z;
            kernel_axpby(N1, p, z, 0, 1, axpby_twcl);
        } else {
            //double beta = po / po_0;
            //kernel_axpby(N1, p, z, beta, 1);
            kernel_axpby(N1, p, z, po / po_0, 1, axpby_twcl);
        }
        kernel_SpMV(A, IA, JA, p, q, N1, SpMV_twcl);
        double alpha = po / kernel_dot(N1, p, q, dot_twcl);
        kernel_axpby(N1, x, p, 1, alpha, axpby_twcl);
        kernel_axpby(N1, r, q, 1, -alpha, axpby_twcl);
        po_0 = po;

		if (debug) {
    		double res_i = kernel_L2_ResidentialRate(A, IA, JA, x, b, N1);
        	std::cout << "After iter " << iterations << " res: " << res_i;
			std::cout << "\n";
		}

        // if (po < tol || iterations >= MAXITER)
		if (po < tol)
            convergence = true;
        else
            iterations += 1;

    } while (!convergence);
    twcl += omp_get_wtime() - tbeg;

    res = kernel_L2_ResidentialRate(A, IA, JA, x, b, N1);

	if (debug) {
    	std::cout << "CG time + res evaluate: " << omp_get_wtime() - tbeg;
    	std::cout << "\nX: ";
    	for (int i = 0; i < N1; i++)
    	    std::cout << " " << x[i];
    	std::cout << "\nres: " << res << "\n";
	}

//	size_t extra_array = 0;
//	for (int i = 0; i < N1; i++) {
//		extra_array += sizeof r[i];
//		extra_array += sizeof z[i];
//		extra_array += sizeof p[i];
//		extra_array += sizeof q[i];
//	}
//	std::cout << "Extra: " << extra_array << "\n";

    delete []r;
    delete []z;
    delete []p;
    delete []q;

    return;
}

// By already given IA and JA generate (fill) matrix A (CSR format)
void Matrices::fill(double &twcl, bool debug) {
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
                sum += fabs(A[_j]);
            }
        }
        A[diag] = 1.234*sum;
        rev_diag[i] = 1 / A[diag]; // M^(-1)
        b[i] = sin(i);
    }
    twcl += omp_get_wtime() - tbeg;

//	std::cout << "Extra: " << sizeof(double) + sizeof(int) << "\n";

	if (debug)
    	std::cout << "Fill time: " << omp_get_wtime() - tbeg << "\n";

    return;
}

