#include "Graph.h"

// simple destructor. delete all vertices of graph
Graph::~Graph()
{
	std::cout << "WOWO\n";
	delete []Vertices;
	std::cout << "WOWO\n";
//    for (int i = 0; i < N; i++) {
//        Vertex *tmp = first;
//        first = first->next;
//        delete tmp;
//    }
}

// Save in matrix structure IA and JA for this graph
void Graph::get_matrices(Matrices *obj, double &twcl, int trn, bool debug) {
	// devide: every thread work with = N/trn + N%trn > 0 ? 1: 0 lines
	if (N <= trn) {
		for (int i = 0; i < N; i++) {
			//Vertex *cur = get_vertex(i);
			Vertex cur = Vertices[i];
			obj->IA[i+1] = obj->IA[i] + cur.count;
			//Node *connect = cur->neighbours;
	
//			int j = 0;
//	        while(connect) {
//	            //Vertex *tmp = connect->vertex;
//				obj->JA[obj->IA[i] + j] = connect->index;
//				++j;
//				connect = connect->next;
//			}
			for (int j = 0; j < cur.count; j++) {
				obj->JA[obj->IA[i] + j] = cur.neighbours[j];
			}
	
	    }
	    return;
	}
	
	//positions
	int thread_lines = N/trn + ((N%trn > 0)?1:0);
	int *thread_start = new int[trn+1];
	for (int i = 0; i < trn; i++) {
		thread_start[i] = i*thread_lines;
	}
	thread_start[trn] = N;

	// shift related thread
	int *shift = new int[trn + 1];		

	shift[0] = 0;
	obj->IA[0] = 0; 
	obj->JA[0] = 0;
	// many level decomposition
	#pragma omp parallel for
	// evaluate inside every thread like if it was not related matrices
	for (int i = 0; i < trn; i++) {
//		Vertex *cur = get_vertex(thread_start[i]);
//		Vertex cur = Vertices[thread_start[i]];
//		std::cout << "\nQQ " << thread_start[i] << "\n";
//		std::cout << "\nI " << Vertices[thread_start[i]].count;
		obj->IA[thread_start[i]+1] = Vertices[thread_start[i]].count;
		shift[i+1] = Vertices[thread_start[i]].count;

		for (int ie = thread_start[i] + 1; ie < thread_start[i+1]; ie++) {
			//Vertex *cur_ie = get_vertex(ie);
		//	Vertex cur_ie = Vertices[thread_start[ie]];
//			std::cout << "\n1QQ " << ie << "\n";
//			std::cout << "\n1I " << Vertices[ie].count;
			obj->IA[ie+1] = obj->IA[ie] + Vertices[ie].count;
			shift[i+1] += Vertices[ie].count;	
		}
	}

//	std::cout << "\nQQQQ\n";
	for (int i = 1; i < trn; i++)
		shift[i+1] += shift[i];
//	std::cout << "\nQQQQ\n";

	#pragma omp parallel for
	for (int i = 0; i < trn; i++) {
		for (int ie = thread_start[i]; ie < thread_start[i+1]; ie++) {
			obj->IA[ie + 1] += shift[i];
		}
	}

	delete [] shift;
	delete [] thread_start;

//	std::cout << "\nQQQQ\n";
	#pragma omp barrier
	#pragma omp parallel for
	for (int i = 0; i < N; i++) {
	//	Vertex cur = Vertices[i];
	//	Vertex *cur = get_vertex(i);
	//	Node *connect = cur->neighbours;
//		int j = 0;
//		int count = Vertices[i].count;
//		std::cout << "! " << count << "\n";
		for (int j = 0; j < Vertices[i].count; j++) {
//			std::cout << "IIIIIIIIIIIIIIIIIIIIIIII " << Vertices[i].neighbours[j] << "\n";								
//			std::cout << "JA " << obj->IA[i] + j << "\n";
			//int index = obj->IA[i] + j;
			obj->JA[obj->IA[i] + j] = Vertices[i].neighbours[j];
//			std::cout << "RESULT JA " << obj->JA[index] << "\n";
		}
//	    while(connect) {
//	    	//Vertex *tmp = connect->vertex;
//			//obj->JA[obj->IA[i] + j] = tmp->index;
//			obj->JA[obj->IA[i] + j] = connect->index;
//			++j;
//			connect = connect->next;
//		}
	}
//	obj->E1 = obj->IA[N];
//	std::cout << "\nQQQQ\n";
//	std::cout << "E! " << obj->E1 << "\n";
//	for (int i = 0; i < obj->E1; i++) {                                     
//        std::cout << obj->JA[i] << " : ";                                 
//	}
	
	//delete [] thread_start;
	//delete [] shift;

//	std::cout << "\n5QQQQ\n";

	return;

//    int *pos_i = new int[N+1];
//    pos_i[0] = 0;
//    int line_pos[N];
//    int line[N][N];
//    double tbeg = omp_get_wtime();
//    #pragma omp parallel for
//    for (int i = 0; i < N; i++) {
//        Vertex *cur = get_vertex(i);
//        Node *connect = cur->neighbours;
//        line_pos[i] = 0;
//        while(connect) {
//            Vertex *tmp = connect->vertex;
//            line[i][line_pos[i]] = tmp->index;
//            line_pos[i] ++;
//            connect = connect->next;
//        }
//    }
//
//    for (int i = 1; i < N; i++) {
//        pos_i[i] = pos_i[i-1] + line_pos[i-1];
//    }
//
//    #pragma omp parallel for
//    for (int i = 0; i < N; i++) {
//        obj->IA[i] = pos_i[i];
//        for (int j = 0; j < line_pos[i]; j++) {
//            obj->JA[pos_i[i] + j] = line[i][j];
//        }
//    }
//    obj->IA[N] = E;
//
//    twcl += omp_get_wtime() - tbeg;
//	if (debug) {
//    	std::cout << "Get_Matrices time: " << twcl << "\n";
//	}

}

// Only for debug
//void Graph::out_neighbour(Vertex *v) {
//	std::cout << "Out_neighbours for vertex: " << v->index << "\n";
//    Node **list = &(v->neighbours);
//    while (*list) {
//        //Vertex *tmp = (*list)->vertex;
//		//std::cout << "INDEX " << tmp->index << "\n";
//		std::cout << "INDEX " << (*list)->index << "\n";
//        list = &(*list)->next;
//    }
//	std::cout << "End out_neighbours for vertex: " << v->index << "\n";
//    return;
//}

//// add neighbours in ordered list
////void Graph::add_neighbour(Vertex *v, Vertex *neighbour) {
//void Graph::add_neighbour(Vertex *v, int neighbour) {
//    Node **list = &(v->neighbours);
//    while(*list) {
//        int tmp_smaller = (*list)->index;
//        int tmp_bigger = 0;
//        if ((*list)->next)
//            tmp_bigger = ((*list)->next)->index;
//
//		// neighbour vertex already added
//        if (tmp_smaller == neighbour)
//            return;
//		// neighbour vertex already added
//        if (tmp_bigger != 0 && tmp_bigger == neighbour)
//            return;
//
//        if (neighbour < tmp_smaller) {
//            //cout << "INPUT SMALLER THAN EXIST " << neighbour << "\n";
//            Node *tmp = new(struct Node);
//            tmp->index = neighbour;
//            tmp->next = *list;
//            *list = tmp;
//            ++E;
//			++v->count;
//            return;
//        } else if (neighbour > tmp_smaller && tmp_bigger != 0 &&
//					neighbour < tmp_bigger) {
//            //cout << "INPUT THIS " << neighbour << " BETWEEN " << tmp_smaller << " AND " << tmp_bigger << "\n";
//            Node *tmp = new(struct Node);
//            tmp->index = neighbour;
//            tmp->next = (*list)->next;
//            (*list)->next = tmp;
//            ++E;
//			++v->count;
//            return;
//        } else if (tmp_bigger == 0) {
//            //cout << "INPUT BIGGER THAN EXIST " << neighbour << "\n";
//            Node *tmp = new(struct Node);
//            tmp->index = neighbour;
//            tmp->next = NULL;
//            (*list)->next = tmp;
//            ++E;
//			++v->count;
//            return;
//        }
//        list = &(*list)->next;
//        if (*list)
//            list = &(*list)->next;
//    }
//
//	// Add first neighbour in list
//    *list = new(struct Node);
//    (*list)->index = neighbour;
//    (*list)->next = NULL;
//    ++E;
//	++v->count;
//    return;
//}

// return vertex object by index 
// return NULL if vertex with given index doesn't exist
//Vertex* Graph::get_vertex(int index) {
//    Vertex *temp = first;
//    while(temp != NULL) {
//        if (temp->index == index)
//            return temp;
//        temp = temp->next;
//    }
//    return NULL;
//}

// add new vertex in graph
//Vertex* Graph::add_vertex(int index) {
//void Graph::add_vertex(int index) {
    //Vertex *item = get_vertex(index);
	//Vertex * item;
   // if (item == NULL) {
//    item = new(struct Vertex);
//    item->index = index;
////    item->value = 1;
//	item->count = 0;
//    item->neighbours = NULL;
//    item->next = NULL;
//
//    Vertex **current = &first;
//    while(*current) {
//        current = &(*current)->next;
//    }
//    *current = item;
//    //}
//    return item;
//}

//// Add this  __
////          |__|
//void Graph::add_quadrangle(int pos, int Nx, bool debug) {
//    // add all verteces     
//    Vertex *LeftTop     = add_vertex(pos);
//    Vertex *RightTop    = add_vertex(pos + 1);
//    Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
//    Vertex *RightBottom = add_vertex(pos + Nx + 2);
//
//    // add adjacent
//    add_neighbour(LeftTop, LeftTop);
//    add_neighbour(LeftTop, RightTop);
//    add_neighbour(LeftTop, LeftBottom);
//
//    add_neighbour(RightTop, RightTop);
//    add_neighbour(RightTop, LeftTop);
//    add_neighbour(RightTop, RightBottom);
//
//    add_neighbour(LeftBottom, LeftBottom);
//    add_neighbour(LeftBottom, LeftTop);
//    add_neighbour(LeftBottom, RightBottom);
//
//    add_neighbour(RightBottom, RightBottom);
//    add_neighbour(RightBottom, LeftBottom);
//    add_neighbour(RightBottom, RightTop);
//	
//	debug = true;
//    if (debug) {
//		out_neighbour(LeftTop);
//    	out_neighbour(RightTop);
//    	out_neighbour(LeftBottom);
//    	out_neighbour(RightBottom);
//   	}
//    return;
//}

//// Add this: __
////          | /|
////          |/_|
//void Graph::add_2triangles(int pos, int Nx, bool debug) {
//    // add all verteces     
//    Vertex *LeftTop     = add_vertex(pos);
//    Vertex *RightTop    = add_vertex(pos + 1);
//    Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
//    Vertex *RightBottom = add_vertex(pos + Nx + 2);
//
//    // add adjacent
//    add_neighbour(LeftTop, LeftTop);
//    add_neighbour(LeftTop, RightTop);
//    add_neighbour(LeftTop, LeftBottom);
//
//    add_neighbour(RightTop, RightTop);
//    add_neighbour(RightTop, LeftBottom);
//    add_neighbour(RightTop, LeftTop);
//    add_neighbour(RightTop, RightBottom);
//
//    add_neighbour(LeftBottom, LeftBottom);
//    add_neighbour(LeftBottom, LeftTop);
//    add_neighbour(LeftBottom, RightBottom);
//    add_neighbour(LeftBottom, RightTop);
//
//    add_neighbour(RightBottom, RightBottom);
//    add_neighbour(RightBottom, LeftBottom);
//    add_neighbour(RightBottom, RightTop);
//
//	debug = true;
//    if (debug) {
//		out_neighbour(LeftTop);
//    	out_neighbour(RightTop);
//    	out_neighbour(LeftBottom);
//    	out_neighbour(RightBottom);
//   	}
//	return;
//}

void Graph::create_graph(int Nx, int Ny, int K1, int K2, double &twcl, \
						bool debug) {
	double tbeg = omp_get_wtime();
	N = (Nx + 1)*(Ny + 1);
//	for (int i = 0; i < N; i++) {
//		std::cout << "Vertices: " << i << " count: " << Vertices[i].count << "\n";	
//	}
	#pragma omp parallel for
	for (int i = 0; i < (Ny + 1); i++) {
	    for (int j = 0; j < (Nx + 1); j++) {
	        int pos = i*(Nx+1) + j; 
//			Vertex *vertex = add_vertex(pos);
// 			Vertices[pos].neighbours[Vertices[pos].count++] = index;	

			Vertices[pos].count = 0;
			//left top corner
			if (pos == 0) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
//				add_neighbour(vertex, pos + 1);
				//add_neighbour(vertex, pos + Nx + 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}

			int top_cell = i * Nx + j - 1;
			int down_cell = top_cell - Nx + 1;

			//right top corner
			if (pos == Nx) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
//				add_neighbour(vertex, pos - 1);
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
				//	add_neighbour(vertex, pos + Nx);	
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
//				add_neighbour(vertex, pos + Nx + 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}
			//left bottom corner
			if (pos == (Nx+1)*Ny) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
//				add_neighbour(vertex, pos - Nx - 1);
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
					//add_neighbour(vertex, pos - Nx);	
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
//				add_neighbour(vertex, pos + 1);
				continue;
			}
			//right bottom corner
			if (pos == N - 1) {
				//add_neighbour(vertex, pos - Nx - 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
//				add_neighbour(vertex, pos - 1);
				continue;
			}

			// top edge (not corner)
			if (i == 0 && j != 0 && j != Nx) {
				//add_neighbour(vertex, pos - 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
				//add_neighbour(vertex, pos + 1);
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
					//add_neighbour(vertex, pos + Nx);	
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
				//add_neighbour(vertex, pos + Nx + 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}

			// left edge (not corner)
			if (j == 0 && i != 0 && i != Ny) {
				//add_neighbour(vertex, pos - Nx - 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
					//add_neighbour(vertex, pos - Nx);	
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
//				add_neighbour(vertex, pos + 1);
//				add_neighbour(vertex, pos + Nx + 1);
				continue;
			}
			
			// right edge (not corner)
			if (j == Nx && i != 0 && i != Ny) {
				//add_neighbour(vertex, pos - Nx - 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
			//	add_neighbour(vertex, pos - 1);
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
					//add_neighbour(vertex, pos + Nx);	
				}
			//	add_neighbour(vertex, pos + Nx + 1);
				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;	
			}

			// bottom edge (not corner)
			if (i == Ny && j != 0 && j != Nx) {
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				//add_neighbour(vertex, pos - Nx - 1);
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
					//add_neighbour(vertex, pos - Nx);	
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
//				add_neighbour(vertex, pos - 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
//				add_neighbour(vertex, pos + 1);
				continue;
			}

			// center
			if (i != 0 && i != Ny && j != 0 && j != Nx) {
				//add_neighbour(vertex, pos - Nx - 1);
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx - 1;	
				if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
				//	add_neighbour(vertex, pos - Nx);	
 					Vertices[pos].neighbours[Vertices[pos].count++] = pos - Nx;	
				}
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos - 1;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos;	
 				Vertices[pos].neighbours[Vertices[pos].count++] = pos + 1;	
				//add_neighbour(vertex, pos - 1);
				//add_neighbour(vertex, pos + 1);
				if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
				//	add_neighbour(vertex, pos + Nx);	
					Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx;	
				}
			//	add_neighbour(vertex, pos + Nx + 1);
				Vertices[pos].neighbours[Vertices[pos].count++] = pos + Nx + 1;	
				continue;
			}

			// for A2
		//	int top_cell = i * Nx + j;
		//	// top diag
		//	if ((top_cell % (K1 + K2)) >= K1 && j != Nx && i != Ny) {
		//		add_neighbour(vertex, pos + Nx + 2);	
		//	}

		//	int down_cell = top_cell - Nx - 1;
		//	// down diag
		//	if (j != 0 && i != 0 && (down_cell % (K1 + K2)) >= K1) {
		//		add_neighbour(vertex, pos - Nx - 2);	
		//	}
			// top diag
	//		int top_cell = i * Nx + j - 1;
	//		if ((top_cell % (K1 + K2)) >= K1 && j != 0 && i != Ny) {
	//			add_neighbour(vertex, pos + Nx);	
	//		}

	//		// down diag
	//		int down_cell = top_cell - Nx + 1;
	//		if (j != Nx && i != 0 && (down_cell % (K1 + K2)) >= K1) {
	//			add_neighbour(vertex, pos - Nx);	
	//		}

	//		add_neighbour(vertex, pos);
			//out_neighbour(vertex);
		}
	}

	//evaluat E
	//count of K1 + k2 repeats
//	int div = (Nx * Ny) / (K1 + K2);
//	int mod = (Nx * Ny) % (K1 + K2);
//	E = div *2 * K2 + 8 + (Nx-1)*3*2 + (Ny-1)*3*2 + (Nx-1)*(Ny-1)*4;
//	if (mod > K1)
//		E = E + (mod - K1)*2;
	E = 0;
	#pragma omp parallel for reduction(+:E)
	for (int i = 0; i < N; i++) {
		E += Vertices[i].count;
	}

	std::cout << "\n\nMEOW " << E << "\n";
//	#pragma omp parallel for
//	for (int i = 0; i < Ny; i++) {
//	    for (int j = 0; j < Nx; j++) {
//	        int pos = i*(Nx+1) + j; 
//	        if (((pos)%(K1+K2+1)) < K1) {
//	            #pragma omp critical
//	            {                   
//	                add_quadrangle(pos, Nx);
//	            }
//	        } else {
//	            #pragma omp critical
//	            {                   
//	                add_2triangles(pos, Nx);
//	            }
//	        }
//	    }
//	}
	twcl += omp_get_wtime() - tbeg;
	if (debug)
		std::cout << "Create Graph time: " << twcl << "\n";

	return;
}
