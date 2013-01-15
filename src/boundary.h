#ifndef BOUNDARY_H
#define BOUNDARY_H

typedef struct s_BOUNDARY * BOUNDARY;

#include "condition.h"
#include "face.h"

#define BOUNDARY_SUCCESS 1
#define BOUNDARY_FAIL 0
#define BOUNDARY_READ_ERROR -1
#define BOUNDARY_MEMORY_ERROR -2
#define BOUNDARY_LOGIC_ERROR -3

#define BOUNDARY_MAX_N_FACES 4

BOUNDARY boundary_new(int index);
int boundary_index(BOUNDARY boundary);
int boundary_read_name(FILE *file, BOUNDARY boundary);
int boundary_read_faces(FILE *file, FACE *face, BOUNDARY boundary);
int boundary_read_condition(FILE *file, BOUNDARY boundary);
CONDITION boundary_condition(BOUNDARY boundary);
void boundary_print(BOUNDARY boundary);
void boundary_free(BOUNDARY boundary);

#endif
