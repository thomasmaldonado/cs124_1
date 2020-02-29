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

// swap two integers with addresses given by ptr1, ptr2
void swap(int* ptr1, int* ptr2) {
	int temp = *ptr1;
	*ptr1 = *ptr2;
	*ptr2 = temp;
	return;
}

// helper function to create an array of consecutive integers
// use for creating labels
int* range(int max) {
	int* range = (int*) malloc(sizeof(int) * max);
	for (int i = 0; i < max; i++)
	{
		range[i] = i;
	}
	return range;
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
		swap(&a, &b);
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


/*
define type Heap to include size of heap, breadth of each level, array of
indices and values, and label/value pairs for constant lookup time

_pos: location in heap
labels[_pos]: label of the vertex at specified location in heap
	label locations are changing a lot
vals[label]: value associated with the label of the vertex at specified location in heap
	values for a given label don't move
*/

typedef struct Heap {
    int size;
    int breadth;
    int* labels; // labels[i] := label of the i^th element in the heap
	float* vals; // values[j] := the value associated with the label j
	int* pos; // position[k] := the position in the heap of the label k
} Heap;

void print_heap(Heap* H){
	printf("Size: %i\n", H -> size);
	printf("Breadth: %i\n", H -> breadth);
	printf("Labels: [ ");
	for(int i = 0; i < H -> size; i++)
	{
		printf("%i ", H -> labels[i]);
	}
	printf("]\n");
	printf("Values: [ ");
	for(int i = 0; i < H -> size; i++)
	{
		printf("%f ", H -> vals[i]);
	}
	printf("]\n");
	printf("Positions: [ ");
	for(int i = 0; i < H -> size; i++)
	{
		printf("%i ", H -> pos[i]);
	}
	printf("]\n\n");
	
}

// return parent_pos of a child in heap
int get_parent_pos(Heap* H, int child_pos) {
	if (child_pos == 0)
	{
		return -1;
	}
	else
	{
		return child_pos / H -> breadth;
	}
}

// return index of specified child given index of its parent 
// (return -1 if child does not exist)
int get_child_pos(Heap* H, int parent_pos, int relative_child_pos) {
    int child_pos = parent_pos * H -> breadth + relative_child_pos + 1;
    if (H -> size <= child_pos)
    {
        return -1;
    }
    else
    {
        return child_pos;
    }
}

void heap_swap(Heap* H, int pos_1, int pos_2) {
	int label_1 = H -> labels[pos_1];
	int label_2 = H -> labels[pos_2];
	swap(&H -> labels[pos_1], &H -> labels[pos_2]);
	return;
}

// reconstruct heap
void min_heapify(Heap* H, int parent_pos) {
	int min_pos = parent_pos;
	int child_pos;
	for (int i = 0; i < H -> breadth; i++)
	{
		child_pos = get_child_pos(H, parent_pos, i);
        if (child_pos == -1)
        {
            break;
        }
		if (H -> vals[H -> labels[child_pos]] < H -> vals[H -> labels[min_pos]])
		{
			min_pos = child_pos;
		}
	}

	if (min_pos != parent_pos)
	{
		// swap locations of the labels in the heap
		heap_swap(H, parent_pos, min_pos);
        min_heapify(H, min_pos);
	}
	return;
}

// build a heap with the given parameters
Heap build_heap(float* vals, int size, int breadth) {
    Heap H;
    H.size = size;
    H.breadth = breadth;
    H.labels = range(H.size);
    H.vals = vals;
    H.pos = H.labels;
    print_heap(&H);
    // bottom-up reconstruction
    for (int i = H.size - 1; i >= 0; i--) // H.size / 2 perhaps
    {
        min_heapify(&H, i);
    }
    return H;
}

// return the label of the minimum element at the head of the heap
int delete_min(Heap* H) {
	printf("Size: %i\n", H -> size);

    H -> size--;

    int min_pos = 0;
    int min_label = H -> labels[min_pos];
    int last_pos = H -> size;
    int last_label = H -> labels[last_pos];
    
    // move the last element in the heap to the top before re-heapification
    // move the first (min) element to the end of heap for consistency
    H -> labels[min_pos] = last_label;
    H -> pos[min_label] = last_pos;
    H -> labels[last_pos] = min_label;
    H -> pos[last_label] = min_pos;

    min_heapify(H, 0);
    return min_label;
}

// insert a label-value pair into the heap
void insert(Heap* H, int insert_label, float insert_val) {
	// don't do anything if current value is smaller than value you're trying to insert
	if (H -> vals[insert_label] <= insert_val)
	{
		return;
	}
	// update value
	H -> vals[insert_label] = insert_val;
	int insert_pos = H -> pos[insert_label];

    // bottom-up reconstruction from insertion position
    int parent_pos = get_parent_pos(H, insert_pos);
    int parent_label = H -> labels[parent_pos];
    int parent_val = H -> vals[parent_label];
    while (parent_val > insert_val && parent_pos != -1)
    {
    	// swap position and label values while necessary
    	heap_swap(H, parent_pos, insert_pos);
    	insert_pos = H -> pos[insert_label];
    	parent_pos = get_parent_pos(H, insert_pos);
		parent_label = H -> labels[parent_pos];
		parent_val = H -> vals[parent_label];
    }
    return;
}

int is_min_heap(Heap* H)
{
	int parent_pos;
	int parent_label;
	int parent_val;
	int child_pos;
	int child_label;
	int child_val;
	for(int i = 0; i < H -> size; i++)
	{
		parent_label = i;
		parent_pos = H -> pos[parent_label];
		parent_val = H -> vals[parent_label];
		for(int j = 0; j < H -> breadth; j++){
			child_pos = get_child_pos(H, parent_pos, j);
			if (child_pos == -1)
			{
				break;
			}
			child_label = H -> labels[child_pos];
			child_val = H -> vals[child_label];
			if (child_val  < parent_val)
			{
				print_heap(H);
				printf("nope! \n\n\n\n\n\n\n\n");
				return 0;
			}
		}



	}
	return 1;
}

// implementation of Prim's alogorithm
// returns total weight of the MST
float prim(Graph G) {
	// intialize arrays of previous nodes and distances from tree to vertices
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
-1
	// initialize heap
	Heap H = build_heap(G.n, 2, range(G.n), dist)

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

int main(int argc, char** argv) {
	int input = atoi(argv[1]);
	int numpoints = atoi(argv[2]);
	int numtrials = atoi(argv[3]);
	int dimension = atoi(argv[4]);
    int breadth = 2;
    float * vals = (float*) malloc(numpoints * sizeof(float));

	// random seed generator
	int seed = time(NULL);
	srand(seed);

    for (int i = 0; i < numpoints; i++)
    {
    	vals[i] = uniform();
    }

	// generate random graph
	// Graph G = rand_graph(numpoints, dimension);
	// print_graph(G);


	// generate heap
    Heap H = build_heap(vals, numpoints, breadth);
    print_heap(&H); 

    printf("Is min heap? %i \n", is_min_heap(&H));

	return 0;
}
