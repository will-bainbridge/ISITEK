#ifndef MEMORY_H
#define MEMORY_H

#include "expression.h"

int ** allocate_integer_matrix(int **old, int height, int width);
double ** allocate_double_matrix(double **old, int height, int width);
char ** allocate_character_matrix(char **old, int height, int width);
void destroy_matrix(void **matrix);
int *** allocate_integer_tensor(int ***old, int height, int width, int depth);
double *** allocate_double_tensor(double ***old, int height, int width, int depth);
void destroy_tensor(void ***tensor);

struct NODE * allocate_nodes(int n_nodes);
void destroy_nodes(int n_nodes, struct NODE *node);

struct FACE * allocate_faces(int n_faces);
struct NODE ** allocate_face_node(struct FACE *face);
struct ELEMENT ** allocate_face_border(struct FACE *face);
int * allocate_face_n_boundaries(struct FACE *face, int n_variables);
struct BOUNDARY *** allocate_face_boundary(struct FACE *face, int n_variables);
double ** allocate_face_x(struct FACE *face, int n_points);
double * allocate_face_w(struct FACE *face, int n_points);
double **** allocate_face_q(struct FACE *face, int n_variables, int *n_basis, int n_points);
void destroy_faces(int n_faces, struct FACE *face, int n_variables);

struct ELEMENT * allocate_elements(int n_elements);
struct FACE ** allocate_element_face(struct ELEMENT *element);
int * allocate_element_orient(struct ELEMENT *element);
int ** allocate_element_unknown(struct ELEMENT *element, int n_variables, int *n_basis);
double ** allocate_element_x(struct ELEMENT *element, int n_points);
double * allocate_element_w(struct ELEMENT *element, int n_points);
double *** allocate_element_p(struct ELEMENT *element, int n_basis, int n_points);
double *** allocate_element_q(struct ELEMENT *element, int n_basis, int n_points);
void destroy_elements(int n_elements, struct ELEMENT *element, int n_variables);

struct BOUNDARY * allocate_boundaries(int n_boundaries);
struct FACE ** allocate_boundary_face(struct BOUNDARY *boundary);
void destroy_boundaries(int n_boundaries, struct BOUNDARY *boundary);

struct TERM * allocate_terms(int n_terms);
int * allocate_term_variable(struct TERM *term);
int * allocate_term_differential(struct TERM *term);
char * allocate_term_method(struct TERM *term);
EXPRESSION * allocate_term_jacobian(struct TERM *term);
void destroy_terms(int n_terms, struct TERM *term);

#endif
