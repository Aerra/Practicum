// Graph implementation

#include "Matrices.h"

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
    void add_neighbour(Vertex *, Vertex *); // Add for one vertex another that adjancent to first given vertex (neighbour)
    void out_neighbour(Vertex *); // ONLY for debug: print all neighbours for given vertex
    Vertex* get_vertex(int);

    void add_quadrangle(int, int, bool=false);
    void add_2triangles(int, int, bool=false);
    int get_N() { return N; };
    int get_E() { return E; };

    void get_matrices(Matrices *, bool=false); // make IA and JA matrices
	void create_graph(int, int, int, int, bool=false);
};


