// CS 124 Programming Assignment 1
// Kelechi Ukah and Thomas Maldonado

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

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
float get_weight(Graph* G, int u, int v) {
	int dim = G -> dim;
	if (dim == 0)
	{
		int index = get_index(u, v, G -> n);
		return G -> edges[index];
	}
	else
	{
		float* u_point = G -> vertices[u];
		float* v_point = G -> vertices[v];
		float dist = 0.0;
		// implement distance formula
		for (int i = 0; i < dim; i++)
		{
			dist += pow(u_point[i] - v_point[i], 2);

		}
    if(dist > dim){
      printf("Bad distance: %f for dimension %i\n", dist, dim);
    }
		return sqrt(dist);
	}
}

float max_weight(Graph* G){
	int d = G -> dim;
	if (d == 0){
		return 1.0;
	}
	else
	{
		return sqrt((float) d);
	}

}

// print everything relevant to the graph
void print_graph(Graph* G) {
	for (int i = 0; i < G -> n; i++)
	{
		for (int j = i+1; j < G -> n; j++)
		{
			printf("Weight between (%i, %i): %f\n", i, j, get_weight(G, i, j));
		}
	}
}

// ----
// heap implementation
typedef struct Heap {
	/*
	define type Heap to include size of heap, breadth of each level, array of
	indices and values, and label/value pairs for constant lookup time

	_pos: location in heap
	labels[_pos]: label of the vertex at specified location in heap
		label locations are changing a lot
	vals[label]: value associated with the label of the vertex at specified location in heap
		values for a given label don't move
	*/
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

// helper function to swap nodes in a heap
void heap_swap(Heap* H, int pos_1, int pos_2) {
	int label_1 = H -> labels[pos_1];
	int label_2 = H -> labels[pos_2];
	swap(&H -> labels[pos_1], &H -> labels[pos_2]);
  swap(&H -> pos[label_1], &H -> pos[label_2]);
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
    H.pos = range(H.size);
    // bottom-up reconstruction
    for (int i = H.size - 1; i >= 0; i--) // H.size / 2 perhaps
    {
        min_heapify(&H, i);
    }
    return H;
}

// return the label of the minimum element at the head of the heap
int delete_min(Heap* H) {
    H -> size--;

    int min_label = H -> labels[0];
    //int last_pos = H -> size;
    //int last_label = H -> labels[last_pos];
    heap_swap(H, 0, H -> size);
    // move the last element in the heap to the top before re-heapification
    // move the first (min) element to the end of heap for consistency
    //H -> labels[0] = last_label;
    //H -> labels[last_pos] = min_label;
    min_heapify(H, 0);
    return min_label;
}

// insert a label-value pair into the heap
void insert(Heap* H, int insert_label, float insert_val) {
	/*
	this case where the insert value is worse than the current value
	is handled by the call of insert() in prims()*/

  if (H -> vals[insert_label] <= insert_val)

	{
    printf("Warning: value increased, no change made");
		return;
	}

	// update value
	H -> vals[insert_label] = insert_val;
	int insert_pos = H -> pos[insert_label];

  // bottom-up reconstruction from insertion position
  int parent_pos = get_parent_pos(H, insert_pos);
  int parent_label = H -> labels[parent_pos];
  float parent_val = H -> vals[parent_label];
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
				return 0;
			}
		}



	}
	return 1;
}

void print_MST(int* prev, int size) {
	for (int i = 0; i < size; i++)
	{
		printf("prev[%i]: %i\n", i, prev[i]);
	}
}
// implementation of Prim's alogorithm
// returns total weight of the MST
float prim(Graph* G, int breadth) {
	// intialize arrays of previous nodes and distances from tree to vertices
	float* dist = (float*) malloc(sizeof(float) * G -> n);
	int* prev = (int*) malloc(sizeof(int) * G -> n);
	int* in_MST = (int*) malloc(sizeof(int) * G -> n);
	float total_weight = 0;

	// initialize heap

	// initalize distance and visited values for starting vertex and start the heap
	dist[0] = 0;
	prev[0] = 0;
	in_MST[0] = 1;
	float max_dist = max_weight(G);
	// initialize distance and visited values for every vertex
	for (int i = 1; i < G -> n; i++)
	{
		dist[i] = max_dist;
		prev[i] = i;
		in_MST[i] = 0;
	}


	Heap H = build_heap(dist, G -> n, breadth);

	// move down the heap and update distances/add edges to MST accordingly
	int v;
	float weight;
	while (H.size > 0)
	{
		// pop off nearest vertex v and commit it to the MST
		// then mark it as in MST and already in the MST
		v = delete_min(&H);
		in_MST[v] = 1;
		total_weight += H.vals[v];

		// get every vertex v not already in the MST, and
		// check every vertex connected to v
		for (int w = 0; w < G -> n; w++)
		{
			// only look at vertices not already in the MST
			if (! in_MST[w]){
				weight = get_weight(G, v, w);
				if (H.vals[w] > weight)
				{
					prev[w] = v;
					insert(&H, w, weight);
				}
			}
		}
	}
	return total_weight;
}

float mean_mst_weight(int breadth, int numpoints, int numtrials, int dimension) {
  Graph G;
  float mean_weight = 0;
  float mst_weight;
  int seed;
  for(int i = 0; i < numtrials; i++){
	  seed = time(NULL);
	  srand(seed);
    G = rand_graph(numpoints, dimension);
    mst_weight = prim(&G, breadth);
    mean_weight += mst_weight;
  }
  mean_weight = mean_weight / numtrials;
  return mean_weight;
}

int main(int argc, char** argv) {
	int breadth = atoi(argv[1]);
	int numpoints = atoi(argv[2]);
	int numtrials = atoi(argv[3]);
	int dimension = atoi(argv[4]);
  if(breadth == 0){
    breadth = 2;
  }
  if(dimension == 1){
    dimension = 0;
  }
  float avg = mean_mst_weight(breadth, numpoints, numtrials, dimension);
	printf("%f %i %i %i\n", avg, numpoints, numtrials, dimension);
	return avg;
}

