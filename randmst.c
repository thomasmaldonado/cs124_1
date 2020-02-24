// CS 124 Programming Assignment 1
// Kelechi Ukah and Thomas Maldonado

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


// structure for graph in adjancency matrix format
// fully connected graph, so matrix is efficient


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



int main(int argc, char** argv) {

	int input = atoi(argv[1]);
	int numpoints = atoi(argv[2]);
	int numtrials = atoi(argv[3]);
	int dimension = atoi(argv[4]);

	// random seed generator
	int seed = time(NULL);
	srand(seed);

	Graph G = rand_graph(numpoints, dimension);
	print_graph(G);
	return 0;
}
