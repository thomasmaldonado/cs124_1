// CS 124 Programming Assignment 1
// Kelechi Ukah and Thomas Maldonado

#include <stdio.h> 
#include <stdlib.h> 


// structure for graph in adjancency matrix format
// fully connected graph, so matrix is efficient
struct Graph
{
	int** adj;
	int n;
}

struct 

// function to create a new adjacency list for a node 
struct Node* newAdjListNode(int dest) 
{ 
    struct adj_list_node* newNode = 
    (struct AdjListNode*) malloc(sizeof(struct AdjListNode)); 
    newNode->dest = dest; 
    newNode->next = NULL; 
    return newNode; 
} 
  



int main(int argc, char** argv) {
	for (int i = 0; i < argc; i++) {
		printf("Argument #%s: %s\n", i, argv[i]);
	}
   return 0;
}