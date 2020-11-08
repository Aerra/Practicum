#include "Graph.h"

// Save in Matrix structure IA and JA for this graph
void Graph::get_matrices(Matrices *obj, double &twcl, int trn, bool debug) {
	// if count of vertices in graph <= count working threads => no ||
    double tbeg = omp_get_wtime();
	if (N <= trn) {
		for (int i = 0; i < N; i++) {
			Vertex cur = Vertices[i];
			obj->IA[i+1] = obj->IA[i] + cur.count;
			for (int j = 0; j < cur.count; j++) {
				obj->JA[obj->IA[i] + j] = cur.neighbours[j];
			}
	
	    }
    	twcl += omp_get_wtime() - tbeg;
	    return;
	}
	
	// devide: every thread work with = N/trn + N%trn > 0 ? 1: 0 lines
	//positions for thread's work area
	int thread_lines = N/trn + ((N%trn > 0)?1:0);
	int *thread_start = new int[trn+1];
	for (int i = 0; i < trn; i++) {
		thread_start[i] = i*thread_lines;
	}
	// just for understanding where work area ends
	thread_start[trn] = N;

	// shift related thread
	int *shift = new int[trn + 1];		

	shift[0] = 0;
	obj->IA[0] = 0; 
	obj->JA[0] = 0;
	// many level decomposition (from lecture 3)
	#pragma omp parallel for
	// evaluate inside every thread like if it was not related matrices
	for (int i = 0; i < trn; i++) {
		obj->IA[thread_start[i]+1] = Vertices[thread_start[i]].count;
		shift[i+1] = Vertices[thread_start[i]].count;

		for (int ie = thread_start[i] + 1; ie < thread_start[i+1]; ie++) {
			obj->IA[ie+1] = obj->IA[ie] + Vertices[ie].count;
			shift[i+1] += Vertices[ie].count;	
		}
	}

	for (int i = 1; i < trn; i++)
		shift[i+1] += shift[i];

	#pragma omp parallel for
	for (int i = 0; i < trn; i++) {
		for (int ie = thread_start[i]; ie < thread_start[i+1]; ie++) {
			obj->IA[ie + 1] += shift[i];
		}
	}

	delete []shift;
	delete []thread_start;

	// i think that can be race without barrier
	#pragma omp barrier
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < Vertices[i].count; j++) {
			obj->JA[obj->IA[i] + j] = Vertices[i].neighbours[j];
		}
	}

    twcl += omp_get_wtime() - tbeg;
	if (debug)
    	std::cout << "Get_Matrices time: " << omp_get_wtime() - tbeg << "\n";
	return;
}

// Only for debug
void Graph::out_neighbour(int pos) {
	if (pos >= N) { return; }
	std::cout << "vertex: " << pos << " -> (";
	for (int i = 0; i < Vertices[pos].count - 1; i++) {
		std::cout << Vertices[pos].neighbours[i] << ", ";
	}
	std::cout << Vertices[pos].neighbours[Vertices[pos].count - 1] << ")\n";
    return;
}

// Fill the Graph object by using given Nx, Ny, K1 and K2 arguments
void Graph::create_graph(int Nx, int Ny, int K1, int K2, double &twcl, \
						bool debug) {
	double tbeg = omp_get_wtime();
	N = (Nx + 1)*(Ny + 1);
	#pragma omp parallel for
	for (int i = 0; i < (Ny + 1); i++) {
	    for (int j = 0; j < (Nx + 1); j++) {
			// current position	
	        int pos = i*(Nx+1) + j; 
			Vertices[pos].count = 0;

			// fill neighbours in already sorted order (ASC)
			// left top corner
			if (pos == 0) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}

			// for better understanding
			//     +--+
			//     | /|
			//     |/ | - down cell
			//  +--X--+
			//  | /| 
			//  |/ | - top cell
			//  +--+
			//  X - our pos (position)
			//  down because X in the bottom of triangle (/)
			//  top in opposite down
			int top_cell = i * Nx + j - 1;
			int down_cell = top_cell - Nx + 1;

			//right top corner
			if (pos == Nx) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}
			//left bottom corner
			if (pos == (Nx+1)*Ny) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
				continue;
			}
			//right bottom corner
			if (pos == N - 1) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
				continue;
			}

			// top edge (not corner)
			if (i == 0 && j != 0 && j != Nx) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}

			// left edge (not corner)
			if (j == 0 && i != 0 && i != Ny) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}
			
			// right edge (not corner)
			if (j == Nx && i != 0 && i != Ny) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;	
			}

			// bottom edge (not corner)
			if (i == Ny && j != 0 && j != Nx) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
				continue;
			}

			// center
			if (i != 0 && i != Ny && j != 0 && j != Nx) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}
		}
	}

	//evaluate E
	E = 0;
	#pragma omp parallel for reduction(+:E)
	for (int i = 0; i < N; i++) {
		E += Vertices[i].count;
	}

	twcl += omp_get_wtime() - tbeg;
	if (debug) {
		std::cout << "Create Graph time: " << omp_get_wtime() - tbeg << "\n";
		for (int i = 0; i < N; i++) {
			out_neighbour(i);
		}
	}

	return;
}
