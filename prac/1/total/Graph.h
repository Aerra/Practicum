// Graph implementation

#include "Matrices.h"
#include "omp.h"

struct Vertex {
	int count = 0;
	int neighbours[10];
};

class Graph {
    int N; // count of verteces
    int E; // already x2 because current graph is directed graph
    Vertex *Vertices; // array of vertices in graph
public:
    Graph(int N) {
		N = N;
		E = 0;
		Vertices = new Vertex[N];
	};
    ~Graph() { delete []Vertices; };
	void out_neighbour(int); // ONLY for debug: print all neighbours for given vertex
    int get_N() { return N; };
    int get_E() { return E; };

    void get_matrices(Matrices *, double &, int, bool=false); // make IA and JA matrices
	void create_graph(int, int, int, int, double &, bool=false);
};


