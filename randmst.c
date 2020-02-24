// CS 124 Programming Assignment 1
// Kelechi Ukah and Thomas Maldonado

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MAX_INT 2147483647


typedef struct Graph
{
	int n;
	int dim;
	float* edges;
	float** vertices;
} Graph;

// generate uniformly distributed float
float uniform() {
	return (float) rand() / RAND_MAX;
}

// generate a random vertex given an integer number of dimensions
float* rand_vertex(int dim) {
	float* vertex = (float*)malloc(sizeof(float) * dim);
	for (int i = 0; i < dim; i++)
	{
		vertex[i] = uniform();
	}
	return vertex;
}

Graph rand_graph(int n, int dim) {
    if (dim < 0)
    {
        printf("Warning: dimension %i < 0. Reassigned to dimension 0. You should probably change the number of dimensions. I think it's a good idea. What do you think?\n", dim);
        dim = 0;
    }

    Graph G;
    G.n = n;
    G.dim = dim;

    if (dim == 0)
    {
        int m = (n*n - n)/2;
        float* edges = (float*) malloc(m * sizeof(float)); // initialize a space of n_choose_2 edge values

        for (int i = 0; i < m; i++)
        {
            edges[i] = uniform();
        }
        G.edges = edges;
    }
    else
    {
        float** vertices = (float**)malloc(sizeof(float) * dim * n);
        for(int i = 0; i < n; i++)
        {
            vertices[i] = rand_vertex(dim);
        }
        G.vertices = vertices;
    }
    return G;
}

// calculate the position in the list of edges given
int get_index(int a, int b, int n) {
	// swap so that a < b
	// to be consistent in finding location in list of edges
	if (a > b)
	{
		a = a + b;
		b = a - b;
		a = a - b;
	}
	return a*n - (a*(a-1)/2) + b - 2*a - 1;
}

// for dim > 0:
// calculate weight between vertices u and v
float weight(Graph G, int u, int v) {
	int dim = G.dim;
	if (dim == 0)
	{
		int index = get_index(u, v, G.n);
		return G.edges[index];
	}
	else
	{
		float* u_point = G.vertices[u];
		float* v_point = G.vertices[v];
		float dist = 0;
		// implement distance formula
		for (int i = 0; i < dim; i++)
		{
			dist += pow(u_point[i] - v_point[i], 2);
		}
		return sqrt(dist);
	}
}

// print everything relevant to the graph
void print_graph(Graph G) {
	for (int i = 0; i < G.n; i++)
	{
		for (int j = i+1; j < G.n; j++)
		{
			if (G.dim > 0)
			{
				printf("Coordinates of %i:\n( ", i);
				for (int k = 0; k < G.dim; k++)
				{
					printf("%f ", G.vertices[i][k]);
				}
				printf(")\n");
				printf("Coordinates of %i:\n( ", j);
				for (int k = 0; k < G.dim; k++)
				{
					printf("%f ", G.vertices[j][k]);
				}
				printf(")\n");
			}
			printf("Weight between (%i, %i): %f\n\n", i, j, weight(G, i, j));
		}
	}
}

// define type Heap to include size of heap, breadth of each level, array of keys and values, along with array of pointers
// to key locations for constant lookup time
typedef struct Heap {
    int size;
    int breadth;
    int* keys;
	float* vals;
	float** ptrs;
} Heap;

// return index of parent value in heap
int parent(Heap H, int child) {
	return child / H.breadth;
}

// return index of specified child given parent index (return -1 if child does not exist)
int child(Heap H, int parent, int child) {
    int child_index = parent * H.breadth + child;
    if (H.size - 1 < child_index)
    {
        return -1;
    }
    else
    {
        return child_index;
    }
}

// reconstruct heap
void min_heapify(Heap H, int parent) {
	int smallest = parent;
	for (int i = 0; i < H.breadth; i++)
	{
		int child_index = child(H, parent, i);
        if (child_index == -1)
        {
            break;
        }
		if (H.vals[child_index] < H.vals[smallest])
		{
			smallest = child_index;
		}
	}
	if (smallest != parent)
	{
        float temp_val = H.vals[parent];
        H.vals[parent] = H.vals[smallest];
        H.vals[smallest] = temp_val;

        int temp_key = H.keys[parent];
        H.keys[parent] = H.keys[smallest];
        H.keys[smallest] = temp_key;


        H.ptrs[parent] = &H.vals[smallest];
        H.ptrs[smallest] = &H.vals[parent];

        min_heapify(H, smallest);
	}
}

Heap build_heap(int size, int breadth, int* keys, float* vals){
    Heap H;
    H.size = size;
    H.breadth = breadth;
    H.keys = keys;
    H.vals = vals;
    H.ptrs = (float**) malloc(sizeof(int*) * size);
    for (int i = 0; i < size; i++)
    {
    	H.ptrs[keys[i]] = &vals[i];
    }
    for (int i = size-1; i >= 0; i--)
    {
        min_heapify(H, i);
    }
    return H;
}

int delete_min(Heap H)
{
    H.size--;
    int min_key = H.keys[0];
    H.keys[0], H.vals[0] = H.keys[H.size], H.vals[H.size];
    min_heapify(H, 0);
    return min_key;
}

// insert a key-value pair into the heap
//void insert(Heap H, int key, float val) {
	// 
	
//}

// implementation of Prim's alogorithm
// returns total weight of the MST
/*
float prim(Graph G) {
	// array of previous nodes and distances from tree to vertices
	float* dist = (float*) malloc(sizeof(float) * G.n);
	int* prev = (int*) malloc(sizeof(int) * G.n);
	int* visited = (int*) malloc(sizeof(int) * G.n);
	float total_weight = 0;

	// initialize distance and visited values for every vertex
	for (int i = 1; i < G.n; i++)
	{
		dist[i] = MAX_INT;
		prev[i] = i;
		visited[i] = 0;
	}

	// 
	Heap H = build_heap(G.n, 2, , dist)

	// initalize distance and visited values for starting vertex and start the heap
	H.vals[0] = 0;
	prev[0] = 0;
	visited[0] = 1;
	insert(H, 0, 0);

	int i = 0;
	while (H.size > 0)
	{
		// pop off nearest vertex v and commit it to the MST
		// then mark it as already in the tree
		int v = delete_min(H);
		visited[v] = 0;
		total_weight += dist[v];
		// get every vertex not already in the tree, and
		// check every vertex connected to v
		for (int w = 0; w < G.n; w++)
		{
			// only look at vertices not already in the MST
			if (!visited[w])
			{
				weight = weight(G, v, w);
				if (dist[w] > weight)
				{
					dist[w] = weight;
					prev[w] = v;
					insert(H, w, weight);
				}
			}
		}
	}
	return total_weight;
}
*/
int main(int argc, char** argv) {
    /*
	int input = atoi(argv[1]);
	int numpoints = atoi(argv[2]);
	int numtrials = atoi(argv[3]);
	int dimension = atoi(argv[4]);

	// random seed generator
	int seed = time(NULL);
	srand(seed);

	Graph G = rand_graph(numpoints, dimension);
	print_graph(G);
    */

    int size = 16;
    int breadth = 2;
    int* keys = malloc(sizeof(int) * size);
    float* vals = malloc(sizeof(float) * size);

    for (int i = 0; i < size; i++)
    {
        keys[i] = i;
        vals[i] = uniform();
        printf("(Key, Value): (%i, %f)\n", keys[i], vals[i]);
    }

    Heap H = build_heap(size, breadth, keys, vals);

    printf("\n\nPointer check:\n\n");
    for (int i = 0; i < size; i++)
    {
    	printf("(Key, Value): %i, %f\n", keys[i], *H.ptrs[i]);
    }

    for (int i = 0; i < size; i++){
        //printf("\nParent %i: %f\n", i, H.vals[i]);
        for (int j = 0; j < breadth; j++)
        {
            int child_index = child(H, i, j);
            if (child_index == -1)
            {
                break;
            }
            //printf("Child %i: %f\n", j, H.vals[child_index]);
            if (H.vals[i] <= H.vals[child_index])
            {
                printf("Success! :D\n");
            }
            else
            {
                printf("Failure! D:\n");
            }
        }
    }
    for (int i = 0; i < size; i++)
    {
        printf("(Key, Value): (%i, %f)\n", H.keys[i], H.vals[i]);
    }

    printf("Pointer check:\n\n");
    for (int i = 0; i < size; i++)
    {
    	printf("Value %f\n", *H.ptrs[i]);
    }

	return 0;
}
