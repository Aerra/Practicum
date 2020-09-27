#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <fstream>
//using namespace std;

struct Node {
	struct Vertex* vertex;
	Node* next;
};

struct Vertex {
	int index;
	Vertex *next;
	struct Node* neighbours; // list of vertices adjacent to this
};

class Graph {
	int N; // count of vertexs
	int E; // already x2 because current graph is directed graph
	Vertex *first; // first vertex in graph
public:
	Graph() : N(0), E(0), first(0) {}
	~Graph();
	Vertex* add_vertex(int);
	void add_neighbour(Vertex *, Vertex *);
	Vertex* get_vertex(int);

	void add_quadrangle(int, int);
	void add_2triangles(int, int);

	void get_matrices(bool, std::string);
};

Graph::~Graph()
{
	for (int i = 0; i < N; i++) {
		Vertex *tmp = first;
		first = first->next;
		delete tmp;
	}
}

void Graph::get_matrices(bool verbose = false, std::string output_file = "") {
	//cout << "START GET MATRICES\n";
	int IA[N+1];	
	int JA[E];
	int pos = 0;
	for (int i = 0; i < N; i++) {
		//TODO What if vertex didn't exist?
		Vertex *cur = get_vertex(i);
		Node *connect = cur->neighbours;
		int line[N];
		int line_pos = 0;
		while(connect) {
			Vertex *tmp = connect->vertex;
			line[line_pos] = tmp->index;
			++line_pos;
			connect = connect->next;
		}				
		//for (int a = 0; a < line_pos; a++){ cout << "GET MA: " << line[a] << "\n"; }

		//TODO Use better sort or maybe resort neighbours when you add it
		for (int start=0; ; start++) {
			bool done = true;
			for(int j = line_pos-2; j >= start; j--) {
				if(line[j+1] < line[j]) {
					int tmp = line[j];
					line[j] = line[j+1];
					line[j+1] = tmp;
					done = false;
				}					
			}			
			if(done)
				break;
		}
		//for (int a = 0; a < line_pos; a++){ cout << "SORT GET MA: " << line[a] << "\n"; }

		IA[i] = pos;
		for (int j = 0; j < line_pos; j++) {
			JA[pos + j] = line[j];
		}			
		pos += line_pos;
	}
	IA[N] = E;

	if (verbose) {
		std::cout << "IA: ";			
		for (int i = 0; i < N+1; i++) {
				std::cout << IA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "JA: ";			
		for (int i = 0; i < E; i++) {
				std::cout << JA[i] << " ";
		}		
		std::cout << "\n";
		std::cout << "N: " << N << "\n";
	}
	if (output_file != (std::string)"") {
		setlocale(LC_ALL, "rus");
		std::ofstream out_file (output_file.c_str());
		if (!out_file.is_open()) {
			std::cout << "Can't open output file " << output_file << "\n"; 
			return;
		}
		out_file << "IA: ";			
		for (int i = 0; i < N+1; i++) {
				out_file << IA[i] << " ";
		}		
		out_file << "\n";
		out_file << "JA: ";			
		for (int i = 0; i < E; i++) {
				out_file << JA[i] << " ";
		}		
		out_file << "\n";
		out_file << "N: " << N << "\n";

		out_file.close();
	}
	return;
}

void Graph::add_neighbour(Vertex *v, Vertex *neighbour) {
	//cout << "Add_neighbours for vertex " << v->index << " this: " << neighbour->index << "\n";
	Node **list = &(v->neighbours);
	while(*list) {
		Vertex *tmp = (*list)->vertex;
		//cout << "Add_neighbour: Exist index: " << tmp->index << "\n";
		if (tmp->index == neighbour->index) {
				return;
		}
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
	add_neighbour(LeftTop, RightTop);
	add_neighbour(LeftTop, LeftBottom);

	add_neighbour(RightTop, LeftTop);
	add_neighbour(RightTop, RightBottom);

	add_neighbour(LeftBottom, LeftTop);
	add_neighbour(LeftBottom, RightBottom);

	add_neighbour(RightBottom, LeftBottom);
	add_neighbour(RightBottom, RightTop);
	return;
}

void Graph::add_2triangles(int pos, int Nx) {
	// add all verteces		
	Vertex *LeftTop     = add_vertex(pos);
	Vertex *RightTop    = add_vertex(pos+1);
	Vertex *LeftBottom  = add_vertex(pos + Nx + 1);
	Vertex *RightBottom = add_vertex(pos + Nx + 2);
	// add adjacent
	add_neighbour(LeftTop, RightTop);
	add_neighbour(LeftTop, LeftBottom);

	add_neighbour(RightTop, LeftBottom);
	add_neighbour(RightTop, LeftTop);
	add_neighbour(RightTop, RightBottom);

	add_neighbour(LeftBottom, LeftTop);
	add_neighbour(LeftBottom, RightBottom);
	add_neighbour(LeftBottom, RightTop);

	add_neighbour(RightBottom, LeftBottom);
	add_neighbour(RightBottom, RightTop);
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

//int stringtoint(const std::string& s) {
//	int n = s.length();
//	std::cout << "LEN " << n << "\n";
//	std::cout << "S " << s << "\n";
//	char str_array[n+1];
//	std::strcpy(str_array, s.c_str());
//	int out = 0;
//	for (int i = 0; i < n; i--) {
//		std::cout << str_array[n-1] << "\n";			
//		out += std::pow(10, i)*str_array[n-i];
//		return out;
//	}
//	return out;
//}


int main(int argc, char *argv[]) {

	std::string empty("");

	//int port = atoi(port1.c_str());
	bool find = false;
	std::string help_s1 = getCmdOption(argc, argv, "--help", find);
	std::string help_s2 = getCmdOption(argc, argv, "-h", find);
	if (find == true) {
		std::cout << "Usage: " << "\n";
		std::cout << "-h [--help] Print this information and return\n";
		std::cout << "-v [--verbose] if exist then print result to stdout\n";
		std::cout << "--Nx=[param] Count of Nx cell (column)\n";
		std::cout << "--Ny=[param] Count of Ny cell (raw)\n";
		std::cout << "--K1=[param] Count of quadrangle\n";
		std::cout << "--K2=[param] Count of diagonally divided quadrilateral\n";
		std::cout << "--output=[filepath] Path to file with result (if not found then verbose param turn on automatically)\n";
		std::cout << "--input=[filepath] Path to file, that contain Nx, Ny, K1 and K2 params (format [Nx]\\n[Ny]\\n[K1]\\n[K2])\n";
		return 0;
	}

	std::string verbose_s1 = getCmdOption(argc, argv, "--verbose", find);
	std::string verbose_s2 = getCmdOption(argc, argv, "-v", find);
	bool verbose = false;
	if (find == true) {
		verbose = true;
	}

	std::string input_file = getCmdOption(argc, argv, "--input=", find);
	int Nx = 0;
	int Ny = 0;
	int K1 = 0;
	int K2 = 0;
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
		in_file.close();
	}

	std::string output_file = getCmdOption(argc, argv, "--output=", find);
	if (find == false) {
		if (verbose == false) {
			std::cout << "No verbose option and no output file => verbose = true\n";
		}
		verbose = true;
	} else {
		setlocale(LC_ALL, "rus");
		std::ofstream out_file (output_file.c_str());
		if (!out_file) {
			std::cout << "Can't open output file " << output_file << "\n"; 
			return -1;
		}
		out_file.close();
	}
	
	Graph A;

	for (int i = 0; i < Ny; i++) {
		for (int j = 0; j < Nx; j++) {
			int pos = i*(Nx+1) + j;
			if (((pos)%(K1+K2+1)) < K1) {
				//cout << "ADD QUADRAT " << pos << "\n";
				A.add_quadrangle(pos, Nx);
			} else {
				//cout << "ADD TREUG " << pos << "\n";
				A.add_2triangles(pos, Nx);
			}
		}	
	}
	A.get_matrices(verbose, output_file);
	return 0;
}
