////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>

#include "expression.h"
#include "memory.h"
#include "structure.h"

////////////////////////////////////////////////////////////////////////////////

int ** allocate_integer_matrix(int **old, int height, int width)
{
	int i;

	int **new = (int **)realloc(old, height * sizeof(int *));
	if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

	new[0] = (int *)realloc(new[0], height * width * sizeof(int));
	if(new[0] == NULL) return NULL;

	for(i = 1; i < height; i++) new[i] = new[i-1] + width;

	return new;
}

double ** allocate_double_matrix(double **old, int height, int width)
{
	int i;

	double **new = (double **)realloc(old, height * sizeof(double *));
	if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

	new[0] = (double *)realloc(new[0], height * width * sizeof(double));
	if(new[0] == NULL) return NULL;

	for(i = 1; i < height; i++) new[i] = new[i-1] + width;

	return new;
}

char ** allocate_character_matrix(char **old, int height, int width)
{
	int i;

	char **new = (char **)realloc(old, height * sizeof(char *));
	if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

	new[0] = (char *)realloc(new[0], height * width * sizeof(char));
	if(new[0] == NULL) return NULL;

	for(i = 1; i < height; i++) new[i] = new[i-1] + width;

	return new;
}

void destroy_matrix(void **matrix)
{
	if(matrix)
	{
		free(matrix[0]);
	}
	free(matrix);
}

int *** allocate_integer_tensor(int ***old, int height, int width, int depth)
{
        int i, j;

        int ***new = (int ***)realloc(old, height * sizeof(int **));
        if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

        new[0] = (int **)realloc(new[0], height * width * sizeof(int *));
        if(new[0] == NULL) return NULL;
	if(old == NULL) new[0][0] = NULL;

        new[0][0] = (int *)realloc(new[0][0], height * width * depth * sizeof(int));
        if(new[0][0] == NULL) return NULL;

        for(i = 1; i < height; i ++) new[i] = new[i-1] + width;
        for(i = 1; i < height; i ++) new[i][0] = new[i-1][0] + width*depth;
        for(i = 0; i < height; i ++) for(j = 1; j < width; j ++) new[i][j] = new[i][j-1] + depth;

        return new;
}

double *** allocate_double_tensor(double ***old, int height, int width, int depth)
{
        int i, j;

        double ***new = (double ***)realloc(old, height * sizeof(double **));
        if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

        new[0] = (double **)realloc(new[0], height * width * sizeof(double *));
        if(new[0] == NULL) return NULL;
	if(old == NULL) new[0][0] = NULL;

        new[0][0] = (double *)realloc(new[0][0], height * width * depth * sizeof(double));
        if(new[0][0] == NULL) return NULL;

        for(i = 1; i < height; i ++) new[i] = new[i-1] + width;
        for(i = 1; i < height; i ++) new[i][0] = new[i-1][0] + width*depth;
        for(i = 0; i < height; i ++) for(j = 1; j < width; j ++) new[i][j] = new[i][j-1] + depth;

        return new;
}

void destroy_tensor(void ***tensor)
{
        if(tensor)
	{
		if(tensor[0])
		{
			free(tensor[0][0]);
		}
		free(tensor[0]);
	}
        free(tensor);
}

////////////////////////////////////////////////////////////////////////////////

struct NODE * allocate_nodes(int n_nodes)
{
	return (struct NODE *)malloc(n_nodes * sizeof(struct NODE));
}

void destroy_nodes(int n_nodes, struct NODE *node)
{
        free(node);
}

////////////////////////////////////////////////////////////////////////////////

struct FACE * allocate_faces(int n_faces)
{
        struct FACE * new = (struct FACE *)malloc(n_faces * sizeof(struct FACE));
	if(new == NULL) return NULL;

	int i;
	struct FACE z = {0,NULL,0,NULL,NULL,NULL,{0.0,0.0},{0.0,0.0},0.0,NULL,NULL,NULL};
	for(i = 0; i < n_faces; i ++) new[i] = z;

	return new;
}

struct NODE ** allocate_face_node(struct FACE *face)
{
	return (struct NODE **)realloc(face->node, face->n_nodes * sizeof(struct NODE *));
}

struct ELEMENT ** allocate_face_border(struct FACE *face)
{
	return (struct ELEMENT **)realloc(face->border, face->n_borders * sizeof(struct ELEMENT *));
}

int * allocate_face_n_boundaries(struct FACE *face, int n_variables)
{
	return (int *)realloc(face->n_boundaries, n_variables * sizeof(struct ELEMENT *));
}

struct BOUNDARY *** allocate_face_boundary(struct FACE *face, int n_variables)
{
	struct BOUNDARY *** new = (struct BOUNDARY ***)realloc(face->boundary, n_variables * sizeof(struct BOUNDARY **));
	if(new == NULL) return NULL;

	int i;

	if(face->boundary == NULL) for(i = 0; i < n_variables; i ++) new[i] = NULL;

	for(i = 0; i < n_variables; i ++)
	{
		if(face->n_boundaries[i])
		{
			new[i] = (struct BOUNDARY **)realloc(new[i], face->n_boundaries[i] * sizeof(struct BOUNDARY *));
			if(new[i] == NULL) return NULL;
		}
		else
		{
			free(new[i]);
			new[i] = NULL;
		}
	}

	return new;
}

double ** allocate_face_x(struct FACE *face, int n_points)
{
	return allocate_double_matrix(face->X, 2, n_points);
}

double * allocate_face_w(struct FACE *face, int n_points)
{
	return (double *)realloc(face->W, n_points * sizeof(double));
}

double **** allocate_face_q(struct FACE *face, int n_variables, int *n_basis, int n_points)
{
	double ****new = (double ****)realloc(face->Q, n_variables * sizeof(double ***));
	if(new == NULL) return NULL;

	int i;

	if(face->Q == NULL) for(i = 0; i < n_variables; i ++) new[i] = NULL;

	for(i = 0; i < n_variables; i ++)
	{
		new[i] = allocate_double_tensor(new[i], n_basis[i], face->n_borders*n_basis[i] + face->n_boundaries[i]*n_points, n_points);
		if(new[i] == NULL) return NULL;
	}

	return new;
}

void destroy_faces(int n_faces, struct FACE *face, int n_variables)
{
        int i, j;
        for(i = 0; i < n_faces; i ++)
        {
                free(face[i].node);
                free(face[i].border);
		free(face[i].n_boundaries);
		if(face[i].boundary) for(j = 0; j < n_variables; j ++) free(face[i].boundary[j]);
		free(face[i].boundary);
		destroy_matrix((void *)face[i].X);
		free(face[i].W);
		if(face[i].Q) for(j = 0; j < n_variables; j ++) destroy_tensor((void *)face[i].Q[j]);
		free(face[i].Q);
        }
        free(face);
}

//////////////////////////////////////////////////////////////////////////////////

struct ELEMENT * allocate_elements(int n_elements)
{
        struct ELEMENT * new = (struct ELEMENT *)malloc(n_elements * sizeof(struct ELEMENT));
	if(new == NULL) return NULL;

	int i;
	struct ELEMENT z = {0,NULL,{0.0,0.0},0.0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
	for(i = 0; i < n_elements; i ++) new[i] = z;

	return new;
}

struct FACE ** allocate_element_face(struct ELEMENT *element)
{
        return (struct FACE **)realloc(element->face, element->n_faces * sizeof(struct FACE *));
}

int ** allocate_element_unknown(struct ELEMENT *element, int n_variables, int *n_basis)
{
	int ** new = (int **)realloc(element->unknown, n_variables * sizeof(int *));
	if(new == NULL) return NULL;

	int i;

	if(element->unknown == NULL) for(i = 0; i < n_variables; i ++) new[i] = NULL;

	for(i = 0; i < n_variables; i ++)
	{
		new[i] = (int *)realloc(new[i], n_basis[i] * sizeof(int));
		if(new[i] == NULL) return NULL;
	}

	return new;
}

double ** allocate_element_x(struct ELEMENT *element, int n_points)
{
	return allocate_double_matrix(element->X, 2, n_points);
}

double * allocate_element_w(struct ELEMENT *element, int n_points)
{
	return (double *)realloc(element->W, n_points * sizeof(double));
}

double *** allocate_element_p(struct ELEMENT *element, int n_basis, int n_points)
{
	return allocate_double_tensor(element->P, n_basis, n_basis, n_points);
}

double *** allocate_element_q(struct ELEMENT *element, int n_basis, int n_points)
{
	return allocate_double_tensor(element->Q, element->n_faces, n_basis, n_points);
}

double *** allocate_element_i(struct ELEMENT *element, int n_variables, int *n_basis, int n_points)
{
	double ***new = (double ***)realloc(element->I, n_variables * sizeof(double **));
	if(new == NULL) return NULL;

	int i;

	if(element->I == NULL) for(i = 0; i < n_variables; i ++) new[i] = NULL;

	for(i = 0; i < n_variables; i ++)
	{
		new[i] = allocate_double_matrix(new[i], n_points, n_basis[i]);
		if(new[i] == NULL) return NULL;
	}

	return new;
}

double ** allocate_element_v(struct ELEMENT *element, int n_basis)
{
	return allocate_double_matrix(element->V, n_basis, element->n_faces);
}

double *** allocate_element_l(struct ELEMENT *element, int n_variables, int *n_basis)
{
	double ***new = (double ***)realloc(element->L, n_variables * sizeof(double **));
	if(new == NULL) return NULL;

	int i;

	if(element->L == NULL) for(i = 0; i < n_variables; i ++) new[i] = NULL;

	for(i = 0; i < n_variables; i ++)
	{
		new[i] = allocate_double_matrix(new[i], n_basis[i], n_basis[i]);
		if(new[i] == NULL) return NULL;
	}

	return new;
}

void destroy_elements(int n_elements, struct ELEMENT *element, int n_variables)
{
        int i, j;
        for(i = 0; i < n_elements; i ++)
        {
                free(element[i].face);
		if(element[i].unknown) for(j = 0; j < n_variables; j ++) free(element[i].unknown[j]);
		free(element[i].unknown);
		destroy_matrix((void *)element[i].X);
		free(element[i].W);
		destroy_tensor((void *)element[i].P);
		destroy_tensor((void *)element[i].Q);
		if(element[i].I) for(j = 0; j < n_variables; j ++) destroy_matrix((void *)element[i].I[j]);
		free(element[i].I);
		destroy_matrix((void *)element[i].V);
		if(element[i].L) for(j = 0; j < n_variables; j ++) destroy_matrix((void *)element[i].L[j]);
		free(element[i].L);
        }
        free(element);
}

//////////////////////////////////////////////////////////////////////////////////

struct BOUNDARY * allocate_boundaries(int n_boundaries)
{
        struct BOUNDARY * new = (struct BOUNDARY *)malloc(n_boundaries * sizeof(struct BOUNDARY));
	if(new == NULL) return NULL;

	int i;
	struct BOUNDARY z = {0,NULL,0,{0,0},0.0};
	for(i = 0; i < n_boundaries; i ++) new[i] = z;

	return new;
}

struct FACE ** allocate_boundary_face(struct BOUNDARY *boundary)
{
	return (struct FACE **)realloc(boundary->face, boundary->n_faces * sizeof(struct FACE *));
}

void destroy_boundaries(int n_boundaries, struct BOUNDARY *boundary)
{
	int i;
	for(i = 0; i < n_boundaries; i ++)
	{
		free(boundary[i].face);
	}
        free(boundary);
}

//////////////////////////////////////////////////////////////////////////////////

struct TERM * allocate_terms(int n_terms)
{
        struct TERM * new = (struct TERM *)malloc(n_terms * sizeof(struct TERM));
	if(new == NULL) return NULL;

	int i;
	struct TERM z = {0,'\0',0.0,0,NULL,NULL,NULL,NULL,NULL,NULL};
	for(i = 0; i < n_terms; i ++) new[i] = z;

	return new;
}

int * allocate_term_variable(struct TERM *term)
{
	return (int *)realloc(term->variable, term->n_variables * sizeof(int));
}

int * allocate_term_differential(struct TERM *term)
{
	return (int *)realloc(term->differential, term->n_variables * sizeof(int));
}

char * allocate_term_method(struct TERM *term)
{
	return (char *)realloc(term->method, (term->n_variables + 1) * sizeof(char));
}

EXPRESSION * allocate_term_weight(struct TERM *term)
{
	EXPRESSION *new = (EXPRESSION *)realloc(term->weight, term->n_variables * sizeof(EXPRESSION));
	if(new == NULL) return NULL;

	int i;
	if(term->weight == NULL) for(i = 0; i < term->n_variables; i ++) new[i] = NULL;

	return new;
}

EXPRESSION * allocate_term_jacobian(struct TERM *term)
{
	EXPRESSION *new = (EXPRESSION *)realloc(term->jacobian, term->n_variables * sizeof(EXPRESSION));
	if(new == NULL) return NULL;

	int i;
	if(term->jacobian == NULL) for(i = 0; i < term->n_variables; i ++) new[i] = NULL;

	return new;
}

void destroy_terms(int n_terms, struct TERM *term)
{
	int i, j;
	for(i = 0; i < n_terms; i ++)
	{
		free(term[i].variable);
		free(term[i].differential);
		free(term[i].method);
		if(term[i].residual) expression_destroy(term[i].residual);
		for(j = 0; j < term[i].n_variables; j ++) if(term[i].weight[j]) expression_destroy(term[i].weight[j]);
		free(term[i].weight);
		for(j = 0; j < term[i].n_variables; j ++) if(term[i].jacobian[j]) expression_destroy(term[i].jacobian[j]);
		free(term[i].jacobian);
	}
	free(term);
}

//////////////////////////////////////////////////////////////////////////////////

EXPRESSION *allocate_initial(int n_variables)
{
	EXPRESSION *new = (EXPRESSION *)malloc(n_variables * sizeof(EXPRESSION));
	if(new == NULL) return NULL;

	int i;
	for(i = 0; i < n_variables; i ++) new[i] = NULL;

	return new;
}

void destroy_initial(int n_variables, EXPRESSION *initial)
{
	int i;
	if(initial)
	{
		for(i = 0; i < n_variables; i ++) if(initial[i]) expression_destroy(initial[i]);
		free(initial);
	}
}

//////////////////////////////////////////////////////////////////////////////////
