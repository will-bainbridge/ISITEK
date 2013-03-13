#ifndef NODE_H
#define NODE_H

#define NODE_SUCCESS 1
#define NODE_FAIL 0
#define NODE_READ_ERROR -1
#define NODE_MEMORY_ERROR -1

typedef struct s_NODE * NODE;

NODE node_new(int index);
int node_index(NODE node);
int node_read_x(FILE *file, NODE node);
const double * node_x(NODE node);
void node_print(NODE node);
void node_free(NODE node);

#endif
