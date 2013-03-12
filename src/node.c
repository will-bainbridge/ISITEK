////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include "node.h"

////////////////////////////////////////////////////////////////////////////////

struct s_NODE
{
	int index;
        double *x;
};

////////////////////////////////////////////////////////////////////////////////

NODE node_new(int index)
{
	NODE node = (NODE)malloc(sizeof(struct s_NODE));
	if(node == NULL) return NULL;

	node->index = index;

	node->x = NULL;

	return node;
}

int node_index(NODE node)
{
	return node->index;
}

////////////////////////////////////////////////////////////////////////////////

int node_read_x(FILE *file, NODE node)
{
	node->x = (double *)malloc(2*sizeof(double));
	if(node->x == NULL) return NODE_MEMORY_ERROR;

	int i;

	for(i = 0; i < 2; i ++)
	{
		if(fscanf(file,"%lf",&(node->x[i])) != 1)
		{
			return NODE_READ_ERROR;
		}
	}

	return NODE_SUCCESS;
}

void node_x(NODE node, double *x)
{
	int i;
	for(i = 0; i < 2; i ++) x[i] = node->x[i];
}

////////////////////////////////////////////////////////////////////////////////

void node_print(NODE node)
{
	printf("node %i\n",node->index);

	int i;
	if(node->x)
	{
		printf("    node->x\n       ");
		for(i = 0; i < 2; i ++) printf("%g ",node->x[i]);
		printf("\n");
	}
}

////////////////////////////////////////////////////////////////////////////////

void node_free(NODE node)
{
	free(node->x);
	free(node);
}

////////////////////////////////////////////////////////////////////////////////
