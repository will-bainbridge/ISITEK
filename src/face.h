#ifndef FACE_H
#define FACE_H

typedef struct s_FACE * FACE;

#include "boundary.h"
#include "element.h"
#include "node.h"
#include "sparse.h"

#define FACE_SUCCESS 1
#define FACE_FAIL 0
#define FACE_READ_ERROR -1
#define FACE_MEMORY_ERROR -2
#define FACE_LOGIC_ERROR -3

FACE face_new(int index);
int face_index(FACE face);
int face_read_nodes(FILE *file, NODE *node, FACE face);
void face_node(FACE face, NODE * node);
int face_add_border(FACE face, ELEMENT element, int index);
void face_border(FACE face, ELEMENT *border);
void face_calculate_size(FACE face);
double face_size(FACE face);
int face_calculate_centre(FACE face);
void face_centre(FACE face, double *centre);
int face_calculate_normal(FACE face);
void face_normal(FACE face, double *normal);
int face_set_boundary(FACE face, BOUNDARY boundary);
int face_calculate_quadrature(FACE face);
int face_n_quadrature(FACE face);
void face_quadrature_x(FACE face, double **x);
void face_quadrature_w(FACE face, double *w);
int face_add_to_system(FACE face, SPARSE system);
void face_print(FACE face);
void face_plot(FACE face);
void face_free(FACE face);

int face_interpolation_start();
void face_interpolation_end();
int face_interpolation_calculate(FACE face);

#endif
