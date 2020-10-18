#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <fstream>
#include <omp.h>
//using namespace std;

struct Matrices {
	int *IA;
	int *JA;
	double *A;
	double *b;
	int N1;
	int E1;
	Matrices(int N, int E) {
		IA = new int[N+1];
		JA = new int[E];
		A = new double[E];
		b = new double[N];
		N1 = N;
		E1 = E;
	};
	~Matrices() {
		delete[] IA;
		delete[] JA;
		delete[] A;
		delete[] b;
	};
	void fill();
};

struct Node {
	struct Vertex* vertex;
	Node* next;
};

struct Vertex {
	int index;
	int value;
	Vertex *next;
	struct Node* neighbours; // list of vertices adjacent to this
};

class Graph {
	int N; // count of verteces
	int E; // already x2 because current graph is directed graph
	Vertex *first; // first vertex in graph
public:
	Graph() : N(0), E(0), first(0) {}
	~Graph();
	Vertex* add_vertex(int);
	void add_neighbour(Vertex *, Vertex *);
	void out_neighbour(Vertex *);
	Vertex* get_vertex(int);

	void add_quadrangle(int, int);
	void add_2triangles(int, int);
	int get_N() { return N; };
	int get_E() { return E; };

	void get_matrices(Matrices *);
};

Graph::~Graph()
{
	for (int i = 0; i < N; i++) {
		Vertex *tmp = first;
		first = first->next;
		delete tmp;
	}
}

void Matrices::fill() {
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
				sum += abs(A[_j]);
			}
		}
		A[diag] = 1.234*sum;
		b[i] = sin(i);
	}
	double twcl = 0;
	twcl += omp_get_wtime() - tbeg;
	std::cout << "TIME 3 " << twcl << "\n";

	return;
}

void Graph::get_matrices(Matrices *obj) {
	//cout << "START GET MATRICES\n";
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
	double twcl = 0;
	twcl += omp_get_wtime() - tbeg;
	std::cout << "TIME 2 " << twcl << "\n";

	//int pos = 0;
	//for (int i = 0; i < N; i++) {
	//	//TODO What if vertex didn't exist?
	//	Vertex *cur = get_vertex(i);
	//	Node *connect = cur->neighbours;
	//	int line[N];
	//	int line_pos = 0;
	//	while(connect) {
	//		Vertex *tmp = connect->vertex;
	//		line[line_pos] = tmp->index;
	//		++line_pos;
	//		connect = connect->next;
	//	}				
	//	//for (int a = 0; a < line_pos; a++){ cout << "SORT GET MA: " << line[a] << "\n"; }

	//	obj->IA[i] = pos;
	//	#pragma omp parallel for
	//	for (int j = 0; j < line_pos; j++) {
	//		obj->JA[pos + j] = line[j];
	//	}			
	//	pos += line_pos;
	//}
	obj->IA[N] = E;

	return;
}

void Graph::out_neighbour(Vertex *v) {
	//cout << "Out_neighbours for vertex " << v->index << "\n";
	Node **list = &(v->neighbours);
	while (*list) {
		//Vertex *tmp = (*list)->vertex;
		//cout << "INDEX " << tmp->index << "\n";
		list = &(*list)->next;
	}
	//cout << "End out_neigbours for vertex: " << v->index << "\n";
	return;				
}

void Graph::add_neighbour(Vertex *v, Vertex *neighbour) {
	//cout << "Add_neighbours for vertex " << v->index << " this: " << neighbour->index << "\n";

	Node **list = &(v->neighbours);
	while(*list) {
		Vertex *tmp_smaller = (*list)->vertex;
		Vertex *tmp_bigger = NULL;
		if ((*list)->next) {
			tmp_bigger = ((*list)->next)->vertex;
		}

		//cout << "Add_neighbour: Exist smaller index: " << tmp_smaller->index << "\n";
		//if (tmp_bigger != NULL) {
		//	cout << "Add_neighbour: Exist bigger index: " << tmp_bigger->index << "\n";
		//}
		if (tmp_smaller->index == neighbour->index) {
			return;
		}
		if (tmp_bigger!= NULL && tmp_bigger->index == neighbour->index) {
			return;
		}

		if (neighbour->index < tmp_smaller->index) {
			//cout << "INPUT SMALLER THAN EXIST " << neighbour->index << "\n";
			Node *tmp = new(struct Node);
			tmp->vertex = neighbour;
			tmp->next = *list;
			*list = tmp;
			++E;
			return;
		} else if (neighbour->index > tmp_smaller->index && tmp_bigger != NULL && neighbour->index < tmp_bigger->index) {
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

	*list = new(struct Node);
	(*list)->vertex = neighbour;
	(*list)->next = NULL;
	++E;
	//cout << "E " << E << "\n";
	return;
}

Vertex* Graph::get_vertex(int index) {
	Vertex *temp = first;
	while(temp != NULL) {
		if (temp->index == index) {
			//cout << "get_vertex: FIND VERTEX " << temp->index << "\n";
			return temp;
		}
		temp = temp->next;
	}
	//cout << "get_vertex: This vertex doesn't exist\n";
	return NULL;
}

Vertex* Graph::add_vertex(int index) {
	Vertex *item = get_vertex(index);
	if (item == NULL) {
		item = new(struct Vertex);
		item->index = index;
		item->value = 1;
		//cout << "INDEX " << item->index << "\n";
		item->neighbours = NULL;
		item->next = NULL;

		Vertex **current = &first;
		while(*current) {
			current = &(*current)->next;
		}
		*current = item;
		++N;
		//cout << "N " << N << "\n";
	}
	return item;
}

void Graph::add_quadrangle(int pos, int Nx) {
	// add all verteces		
	Vertex *LeftTop     = add_vertex(pos);
	Vertex *RightTop    = add_vertex(pos+1);
	Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
	Vertex *RightBottom = add_vertex(pos + Nx + 2);
	// add adjacent
	add_neighbour(LeftTop, LeftTop);
	add_neighbour(LeftTop, RightTop);
	add_neighbour(LeftTop, LeftBottom);
	//out_neighbour(LeftTop);

	add_neighbour(RightTop, RightTop);
	add_neighbour(RightTop, LeftTop);
	add_neighbour(RightTop, RightBottom);
	//out_neighbour(RightTop);

	add_neighbour(LeftBottom, LeftBottom);
	add_neighbour(LeftBottom, LeftTop);
	add_neighbour(LeftBottom, RightBottom);
	//out_neighbour(LeftBottom);

	add_neighbour(RightBottom, RightBottom);
	add_neighbour(RightBottom, LeftBottom);
	add_neighbour(RightBottom, RightTop);
	//out_neighbour(RightBottom);
	return;
}

void Graph::add_2triangles(int pos, int Nx) {
	// add all verteces		
	Vertex *LeftTop     = add_vertex(pos);
	Vertex *RightTop    = add_vertex(pos+1);
	Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
	Vertex *RightBottom = add_vertex(pos + Nx + 2);
	// add adjacent
	add_neighbour(LeftTop, LeftTop);
	add_neighbour(LeftTop, RightTop);
	add_neighbour(LeftTop, LeftBottom);
	//out_neighbour(LeftTop);

	add_neighbour(RightTop, RightTop);
	add_neighbour(RightTop, LeftBottom);
	add_neighbour(RightTop, LeftTop);
	add_neighbour(RightTop, RightBottom);
	//out_neighbour(RightTop);

	add_neighbour(LeftBottom, LeftBottom);
	add_neighbour(LeftBottom, LeftTop);
	add_neighbour(LeftBottom, RightBottom);
	add_neighbour(LeftBottom, RightTop);
	//out_neighbour(LeftBottom);

	add_neighbour(RightBottom, RightBottom);
	add_neighbour(RightBottom, LeftBottom);
	add_neighbour(RightBottom, RightTop);
	//out_neighbour(RightBottom);
}

void print_info(Matrices *obj, bool verbose = false, std::string output_file = "") {
	if (verbose) {
		std::cout << "IA: ";			
		for (int i = 0; i < obj->N1 + 1; i++) {
				std::cout << obj->IA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "JA: ";			
		for (int i = 0; i < obj->E1; i++) {
				std::cout << obj->JA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "N: " << obj->N1 << "\n";
		std::cout << "A: ";			
		for (int i = 0; i < obj->E1; i++) {
				std::cout << obj->A[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "b: ";			
		for (int i = 0; i < obj->N1; i++) {
				std::cout << obj->b[i] << " ";
		}		
		std::cout << "\n";
	}
	if (output_file != (std::string)"") {
		setlocale(LC_ALL, "rus");
		std::ofstream out_file (output_file.c_str());
		if (!out_file.is_open()) {
			std::cout << "Can't open output file " << output_file << "\n"; 
			return;
		}
		out_file << "IA: ";			
		for (int i = 0; i < obj->N1 + 1; i++) {
				out_file << obj->IA[i] << " ";
		}		
		out_file << "\n";
		out_file << "JA: ";			
		for (int i = 0; i < obj->E1; i++) {
				out_file << obj->JA[i] << " ";
		}		
		out_file << "\n";
		out_file << "N: " << obj->N1 << "\n";

		out_file << "A: ";			
		for (int i = 0; i < obj->E1; i++) {
				out_file << obj->A[i] << " ";
		}		
		out_file << "\n";
		out_file << "b: ";			
		for (int i = 0; i < obj->N1; i++) {
				out_file << obj->b[i] << " ";
		}		
		out_file << "\n";

		out_file.close();
	}
}

// ONLY for parse incoming arguments
std::string getCmdOption(int argc, char* argv[], const std::string& option, bool &find)
{
	std::string cmd;
	find = false;
	for( int i = 0; i < argc; ++i)
	{
		std::string arg = argv[i];
		if(0 == arg.find(option))
		{
			find = true;				
		    //std::size_t found = arg.find_last_of(option);
		    cmd = arg.substr(option.length());
		    return cmd;
		}
	}
	return cmd;
}


int main(int argc, char *argv[]) {

	std::string empty("");

	bool find = false;
	std::string help_s = getCmdOption(argc, argv, "-h", find);
	if (find == true) {
		std::cout << "Usage: " << "\n";
		std::cout << "-h Print this information and return\n";
		std::cout << "-v if exist then print result to stdout\n";
		std::cout << "--T=[param] Count of threads\n";
		std::cout << "--Nx=[param] Count of Nx cell (column)\n";
		std::cout << "--Ny=[param] Count of Ny cell (raw)\n";
		std::cout << "--K1=[param] Count of quadrangle\n";
		std::cout << "--K2=[param] Count of diagonally divided quadrilateral\n";
		std::cout << "--output=[filepath] Path to file with result (if not found then verbose param not turn on automatically)\n";
		std::cout << "--input=[filepath] Path to file, that contain Nx, Ny, K1 and K2 params (format [Nx]\\n[Ny]\\n[K1]\\n[K2])\n";
		return 0;
	}

	std::string verbose_s = getCmdOption(argc, argv, "-v", find);
	bool verbose = false;
	if (find == true) {
		verbose = true;
	}

	std::string input_file = getCmdOption(argc, argv, "--input=", find);
	int Nx = 0;
	int Ny = 0;
	int K1 = 0;
	int K2 = 0;
	int T  = 1;
	if (find == false) {
		std::string Nx_s = getCmdOption(argc, argv, "--Nx=", find);
		if (find == false) {
			std::cout << "Can't find Nx , required!\n";
			return -1;
		}
		Nx = atoi(Nx_s.c_str());
		std::string Ny_s = getCmdOption(argc, argv, "--Ny=", find);
		if (find == false) {
			std::cout << "Can't find Ny, required!\n";
			return -1;
		}
		Ny = atoi(Ny_s.c_str());
		std::string K1_s = getCmdOption(argc, argv, "--K1=", find);
		if (find == false) {
			std::cout << "Can't find K1, required!\n";
			return -1;
		}
		K1 = atoi(K1_s.c_str());
		std::string K2_s = getCmdOption(argc, argv, "--K2=", find);
		if (find == false) {
			std::cout << "Can't find K2, required!\n";
			return -1;
		}
		K2 = atoi(K2_s.c_str());
		std::string T_s = getCmdOption(argc, argv, "--T=", find);
		if (find == false) {
			std::cout << "Can't find T, required!\n";
			return -1;
		}
		T = atoi(T_s.c_str());
	} else {
		// Read input file
		setlocale(LC_ALL, "rus");
		std::ifstream in_file (input_file.c_str());
		if (!in_file) {
			std::cout << "Can't open input file " << input_file << "\n"; 
			return -1;
		}
		std::string Nx_s;
	  	getline(in_file, Nx_s);
		if (Nx_s == empty) { 
			std::cout << "Can't find Nx , required!\n";
			return -1;
		}
		Nx = atoi(Nx_s.c_str());
		std::string Ny_s;
	  	getline(in_file, Ny_s);
		if (Ny_s == empty) { 
			std::cout << "Can't find Ny, required!\n";
			return -1;
		}
		Ny = atoi(Ny_s.c_str());
		std::string K1_s;
	  	getline(in_file, K1_s);
		if (K1_s == empty) { 
			std::cout << "Can't find K1, required!\n";
			return -1;
		}
		K1 = atoi(K1_s.c_str());
		std::string K2_s;
	  	getline(in_file, K2_s);
		if (K2_s == empty) { 
			std::cout << "Can't find K2, required!\n";
			return -1;
		}
		K2 = atoi(K2_s.c_str());
		std::string T_s;
	  	getline(in_file, T_s);
		if (T_s == empty) { 
			std::cout << "Can't find T, required!\n";
			return -1;
		}
		T = atoi(T_s.c_str());
		in_file.close();
	}

	std::string output_file = getCmdOption(argc, argv, "--output=", find);
	//if (find == false) {
	//	if (verbose == false) {
	//		std::cout << "No verbose option and no output file => verbose = false \n";
	//	}
	//	verbose = true;
	//} else {
	if (find == true) {
		setlocale(LC_ALL, "rus");
		std::ofstream out_file (output_file.c_str());
		if (!out_file) {
			std::cout << "Can't open output file " << output_file << "\n"; 
			return -1;
		}
		out_file.close();
	}
	
	omp_set_num_threads(T);
	// GENERATE STAGE
	Graph A;
	
	double tbeg = omp_get_wtime();
	#pragma omp parallel for
	for (int i = 0; i < Ny; i++) {
		for (int j = 0; j < Nx; j++) {
			int pos = i*(Nx+1) + j;
			if (((pos)%(K1+K2+1)) < K1) {
				//cout << "ADD QUADRAT " << pos << "\n";
				#pragma omp critical
				{					
					A.add_quadrangle(pos, Nx);
				}
			} else {
				//cout << "ADD TREUG " << pos << "\n";
				#pragma omp critical
				{					
					A.add_2triangles(pos, Nx);
				}
			}
		}	
	}
	double twcl = 0;
	twcl += omp_get_wtime() - tbeg;
	std::cout << "TIME 1 " << twcl << "\n";

	Matrices B(A.get_N(), A.get_E());
	A.get_matrices(&B);

	// FILL STAGE
	B.fill();

	print_info(&B, verbose, output_file);
	return 0;
}
