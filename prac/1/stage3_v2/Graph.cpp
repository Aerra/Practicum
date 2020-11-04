#include "Graph.h"

// simple destructor. delete all vertices of graph
Graph::~Graph()
{
    for (int i = 0; i < N; i++) {
        Vertex *tmp = first;
        first = first->next;
        delete tmp;
    }
}

// Save in matrix structure IA and JA for this graph
void Graph::get_matrices(Matrices *obj, bool debug) {
    int pos_i[N+1];
    pos_i[0] = 0;
    int line_pos[N];
    int line[N][N];
    double tbeg = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        Vertex *cur = get_vertex(i);
        Node *connect = cur->neighbours;
        line_pos[i] = 0;
        while(connect) {
            Vertex *tmp = connect->vertex;
            line[i][line_pos[i]] = tmp->index;
            line_pos[i] ++;
            connect = connect->next;
        }
    }

    for (int i = 1; i < N; i++) {
        pos_i[i] = pos_i[i-1] + line_pos[i-1];
    }

    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        obj->IA[i] = pos_i[i];
        for (int j = 0; j < line_pos[i]; j++) {
            obj->JA[pos_i[i] + j] = line[i][j];
        }
    }
    obj->IA[N] = E;

    double twcl = 0;
    twcl += omp_get_wtime() - tbeg;
	if (debug) {
    	std::cout << "Get_Matrices time: " << twcl << "\n";
	}

    //int pos = 0;
    //for (int i = 0; i < N; i++) {
    //  Vertex *cur = get_vertex(i);
    //  Node *connect = cur->neighbours;
    //  int line[N];
    //  int line_pos = 0;
    //  while(connect) {
    //      Vertex *tmp = connect->vertex;
    //      line[line_pos] = tmp->index;
    //      ++line_pos;
    //      connect = connect->next;
    //  }               
    //  obj->IA[i] = pos;
    //  for (int j = 0; j < line_pos; j++) {
    //      obj->JA[pos + j] = line[j];
    //  }           
    //  pos += line_pos;
    //}
    return;
}

// Only for debug
void Graph::out_neighbour(Vertex *v) {
	std::cout << "Out_neighbours for vertex: " << v->index << "\n";
    Node **list = &(v->neighbours);
    while (*list) {
        Vertex *tmp = (*list)->vertex;
		std::cout << "INDEX " << tmp->index << "\n";
        list = &(*list)->next;
    }
	std::cout << "End out_neigbours for vertex: " << v->index << "\n";
    return;
}

// add neighbours in ordered list
void Graph::add_neighbour(Vertex *v, Vertex *neighbour) {
    Node **list = &(v->neighbours);
    while(*list) {
        Vertex *tmp_smaller = (*list)->vertex;
        Vertex *tmp_bigger = NULL;
        if ((*list)->next)
            tmp_bigger = ((*list)->next)->vertex;

		// neighbour vertex already added
        if (tmp_smaller->index == neighbour->index)
            return;
		// neighbour vertex already added
        if (tmp_bigger != NULL && tmp_bigger->index == neighbour->index)
            return;

        if (neighbour->index < tmp_smaller->index) {
            //cout << "INPUT SMALLER THAN EXIST " << neighbour->index << "\n";
            Node *tmp = new(struct Node);
            tmp->vertex = neighbour;
            tmp->next = *list;
            *list = tmp;
            ++E;
            return;
        } else if (neighbour->index > tmp_smaller->index && \
				tmp_bigger != NULL && neighbour->index < tmp_bigger->index) {
            //cout << "INPUT THIS " << neighbour->index << " BETWEEN " << tmp_smaller->index << " AND " << tmp_bigger->index << "\n";
            Node *tmp = new(struct Node);
            tmp->vertex = neighbour;
            tmp->next = (*list)->next;
            (*list)->next = tmp;
            ++E;
            return;
        } else if (tmp_bigger == NULL) {
            //cout << "INPUT BIGGER THAN EXIST " << neighbour->index << "\n";
            Node *tmp = new(struct Node);
            tmp->vertex = neighbour;
            tmp->next = NULL;
            (*list)->next = tmp;
            ++E;
            return;
        }
        list = &(*list)->next;
        if (*list)
            list = &(*list)->next;
    }

	// Add first neighbour in list
    *list = new(struct Node);
    (*list)->vertex = neighbour;
    (*list)->next = NULL;
    ++E;
    return;
}

// return vertex object by index 
// return NULL if vertex with given index doesn't exist
Vertex* Graph::get_vertex(int index) {
    Vertex *temp = first;
    while(temp != NULL) {
        if (temp->index == index)
            return temp;
        temp = temp->next;
    }
    return NULL;
}

// add new vertex in graph
Vertex* Graph::add_vertex(int index) {
    Vertex *item = get_vertex(index);
    if (item == NULL) {
        item = new(struct Vertex);
        item->index = index;
        item->value = 1;
        item->neighbours = NULL;
        item->next = NULL;

        Vertex **current = &first;
        while(*current) {
            current = &(*current)->next;
        }
        *current = item;
        ++N;
    }
    return item;
}

// Add this  __
//          |__|
void Graph::add_quadrangle(int pos, int Nx, bool debug) {
    // add all verteces     
    Vertex *LeftTop     = add_vertex(pos);
    Vertex *RightTop    = add_vertex(pos + 1);
    Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
    Vertex *RightBottom = add_vertex(pos + Nx + 2);

    // add adjacent
    add_neighbour(LeftTop, LeftTop);
    add_neighbour(LeftTop, RightTop);
    add_neighbour(LeftTop, LeftBottom);

    add_neighbour(RightTop, RightTop);
    add_neighbour(RightTop, LeftTop);
    add_neighbour(RightTop, RightBottom);

    add_neighbour(LeftBottom, LeftBottom);
    add_neighbour(LeftBottom, LeftTop);
    add_neighbour(LeftBottom, RightBottom);

    add_neighbour(RightBottom, RightBottom);
    add_neighbour(RightBottom, LeftBottom);
    add_neighbour(RightBottom, RightTop);
	
    if (debug) {
		out_neighbour(LeftTop);
    	out_neighbour(RightTop);
    	out_neighbour(LeftBottom);
    	out_neighbour(RightBottom);
   	}
    return;
}

// Add this: __
//          | /|
//          |/_|
void Graph::add_2triangles(int pos, int Nx, bool debug) {
    // add all verteces     
    Vertex *LeftTop     = add_vertex(pos);
    Vertex *RightTop    = add_vertex(pos + 1);
    Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
    Vertex *RightBottom = add_vertex(pos + Nx + 2);

    // add adjacent
    add_neighbour(LeftTop, LeftTop);
    add_neighbour(LeftTop, RightTop);
    add_neighbour(LeftTop, LeftBottom);

    add_neighbour(RightTop, RightTop);
    add_neighbour(RightTop, LeftBottom);
    add_neighbour(RightTop, LeftTop);
    add_neighbour(RightTop, RightBottom);

    add_neighbour(LeftBottom, LeftBottom);
    add_neighbour(LeftBottom, LeftTop);
    add_neighbour(LeftBottom, RightBottom);
    add_neighbour(LeftBottom, RightTop);

    add_neighbour(RightBottom, RightBottom);
    add_neighbour(RightBottom, LeftBottom);
    add_neighbour(RightBottom, RightTop);

    if (debug) {
		out_neighbour(LeftTop);
    	out_neighbour(RightTop);
    	out_neighbour(LeftBottom);
    	out_neighbour(RightBottom);
   	}
	return;
}

void Graph::create_graph(int Nx, int Ny, int K1, int K2, bool debug) {
	double tbeg = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < Ny; i++) {
	    for (int j = 0; j < Nx; j++) {
	        int pos = i*(Nx+1) + j; 
	        if (((pos)%(K1+K2+1)) < K1) {
	            #pragma omp critical
	            {                   
	                add_quadrangle(pos, Nx);
	            }
	        } else {
	            #pragma omp critical
	            {                   
	                add_2triangles(pos, Nx);
	            }
	        }
	    }
	}
	double twcl = 0;
	twcl += omp_get_wtime() - tbeg;
	if (debug)
		std::cout << "Create Graph time: " << twcl << "\n";

	return;
}
