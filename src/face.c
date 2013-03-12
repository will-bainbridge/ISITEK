////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "face.h"
#include "memory.h"
#include "quadrature.h"
#include "solver.h"
#include "utility.h"

////////////////////////////////////////////////////////////////////////////////

#include "face.r"

////////////////////////////////////////////////////////////////////////////////

FACE face_new(int index)
{
	FACE face = (FACE)malloc(sizeof(struct s_FACE));
	if(face == NULL) return NULL;

	face->index = index;

	face->node = NULL;
	
	face->n_borders = 0;
	face->border = NULL;

	face->boundary = NULL;

	face->normal = NULL;
	face->centre = NULL;
	face->size = 0.0;

	face->n_quadrature = 0;
	face->X = NULL;
	face->W = NULL;

	face->system_index = -1;

	face->Q = NULL;

	return face;
}

////////////////////////////////////////////////////////////////////////////////

int face_index(FACE face)
{
	return face->index;
}

////////////////////////////////////////////////////////////////////////////////

int face_read_nodes(FILE *file, NODE *node, FACE face)
{
	face->node = (NODE *)malloc(2 * sizeof(NODE));
	if(face->node == NULL) return FACE_MEMORY_ERROR;

	int i, index;

	for(i = 0; i < 2; i ++)
	{
		if(fscanf(file,"%i",&index) != 1)
		{
			return FACE_READ_ERROR;
		}
		face->node[i] = node[index];
	}

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

void face_node(FACE face, NODE * node)
{
	int i;
	for(i = 0; i < 2; i ++) node[i] = face->node[i];
}

////////////////////////////////////////////////////////////////////////////////

int face_add_border(FACE face, ELEMENT element, int index)
{
	if(face->n_borders == 0 && index == 1)
	{
		NODE temp = face->node[0];
		face->node[0] = face->node[1];
		face->node[1] = temp;
	}
	if(face->n_borders == 1 && index == 0)
	{
		return FACE_LOGIC_ERROR;
	}

	face->n_borders ++;

	face->border = (ELEMENT *)realloc(face->border, face->n_borders * sizeof(ELEMENT));
	if(face->border == NULL) return FACE_MEMORY_ERROR;

	face->border[face->n_borders-1] = element;

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int face_n_borders(FACE face)
{
	return face->n_borders;
}

////////////////////////////////////////////////////////////////////////////////

void face_border(FACE face, ELEMENT *border)
{
	int i;
	for(i = 0; i < face->n_borders; i ++) border[i] = face->border[i];
}

////////////////////////////////////////////////////////////////////////////////

void face_calculate_size(FACE face)
{
	int i;

	double x[2][2];
	for(i = 0; i < 2; i ++) node_x(face->node[i],x[i]);

	double dx;
	face->size = 0;
	for(i = 0; i < 2; i ++)
	{
		dx = x[1][i] - x[0][i];
		face->size += dx*dx;
	}

	face->size = sqrt(face->size);
}

////////////////////////////////////////////////////////////////////////////////

double face_size(FACE face)
{
	return face->size;
}

////////////////////////////////////////////////////////////////////////////////

int face_calculate_centre(FACE face)
{
	face->centre = (double *)malloc(2*sizeof(double));
	if(face->centre == NULL) return FACE_MEMORY_ERROR;

	int i;

	double x[2][2];
	for(i = 0; i < 2; i ++) node_x(face->node[i],x[i]);

	for(i = 0; i < 2; i ++)
	{
		face->centre[i] = 0.5*(x[1][i] + x[0][i]);
	}

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

void face_centre(FACE face, double *centre)
{
	int i;
	for(i = 0; i < 2; i ++) centre[i] = face->centre[i];
}

////////////////////////////////////////////////////////////////////////////////

int face_calculate_normal(FACE face)
{
	face->normal = (double *)malloc(2*sizeof(double));
	if(face->normal == NULL) return FACE_MEMORY_ERROR;

	int i;

	double x[2][2];
	for(i = 0; i < 2; i ++) node_x(face->node[i],x[i]);

	face->normal[0] = (- x[1][1] + x[0][1])/face->size;
	face->normal[1] = (+ x[1][0] - x[0][0])/face->size;

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

void face_normal(FACE face, double *normal)
{
	int i;
	for(i = 0; i < 2; i ++) normal[i] = face->normal[i];
}

////////////////////////////////////////////////////////////////////////////////

int face_set_boundary(FACE face, BOUNDARY boundary)
{
	if(face->boundary != NULL) return FACE_LOGIC_ERROR;
	face->boundary = boundary;
	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int face_calculate_quadrature(FACE face)
{
	int i, j;

	face->n_quadrature = solver_n_gauss();
	face->X = matrix_double_new(face->X, 2, face->n_quadrature);
	if(face->X == NULL) return FACE_MEMORY_ERROR;
	face->W = (double *)malloc(face->n_quadrature * sizeof(double));
	if(face->W == NULL) return FACE_MEMORY_ERROR;

	double x[2][2];
	for(i = 0; i < 2; i ++) node_x(face->node[i],x[i]);

	for(i = 0; i < solver_n_gauss(); i ++)
	{
		for(j = 0; j < 2; j ++)
		{
			face->X[j][i] =
				0.5 * (1.0 - quadrature_gauss_location(solver_n_gauss()-1,i)) * x[0][j] +
				0.5 * (1.0 + quadrature_gauss_location(solver_n_gauss()-1,i)) * x[1][j];
		}
		face->W[i] = 0.5 * face->size * quadrature_gauss_weight(solver_n_gauss()-1,i);
	}

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int face_n_quadrature(FACE face)
{
	return face->n_quadrature;
}

////////////////////////////////////////////////////////////////////////////////

void face_quadrature_x(FACE face, double **x)
{
	int i, j;
	for(i = 0; i < face->n_quadrature; i ++)
		for(j = 0; j < 2; j ++)
			x[j][i] = face->X[j][i];
}

////////////////////////////////////////////////////////////////////////////////

void face_quadrature_w(FACE face, double *w)
{
	int i;
	for(i = 0; i < face->n_quadrature; i ++)
		w[i] = face->W[i];
}

////////////////////////////////////////////////////////////////////////////////

int face_add_to_system(FACE face, SPARSE system)
{
	int i, j;

	int n_variables = solver_n_variables(), *n_bases = (int *)malloc(n_variables * sizeof(int)), sum_n_bases = solver_variable_sum_n_bases();
	if(n_bases == NULL) return ELEMENT_MEMORY_ERROR;
	solver_variable_n_bases(n_bases);
	
	int ***unknown;
	unknown = (int ***)malloc(face->n_borders * sizeof(int **));
	if(unknown == NULL) return FACE_MEMORY_ERROR;
	unknown[0] = (int **)malloc(face->n_borders * n_variables * sizeof(int *));
	if(unknown[0] == NULL) return FACE_MEMORY_ERROR;
	unknown[0][0] = (int *)malloc(face->n_borders * sum_n_bases * sizeof(int));
	if(unknown[0][0] == NULL) return FACE_MEMORY_ERROR;
	for(i = 1; i < face->n_borders; i ++) unknown[i] = unknown[i-1] + n_variables;
	for(i = 1; i < face->n_borders; i ++) unknown[i][0] = unknown[i-1][0] + sum_n_bases;
	for(i = 0; i < face->n_borders; i ++) for(j = 1; j < n_variables; j ++) unknown[i][j] = unknown[i][j-1] + n_bases[j-1];

	for(i = 0; i < face->n_borders; i ++) element_unknown(face->border[i],unknown[i]);

	face->system_index = sparse_insert_sub_matrix(system,face->n_borders*sum_n_bases,unknown[0][0],face->n_borders*sum_n_bases,unknown[0][0]);

	free(n_bases);
	free(unknown[0][0]);
	free(unknown[0]);
	free(unknown);

	return face->system_index >= 0 ? FACE_SUCCESS : FACE_LOGIC_ERROR;
}

////////////////////////////////////////////////////////////////////////////////

void face_print(FACE face)
{
	printf("face %i\n",face->index);

	int i, j, k;

	int *n_bases = (int *)malloc(solver_n_variables() * sizeof(int));
	int *interpolation_variable = (int *)malloc(solver_n_interpolations() * sizeof(int));
	int *n_constraints = (int *)calloc(solver_n_variables(), sizeof(int));
	int *constraint_temporary = (int *)malloc(condition_max_n_variables() * sizeof(int));

	solver_variable_n_bases(n_bases);
	solver_interpolation_variable(interpolation_variable);
	if(face->boundary) condition_variable(boundary_condition(face->boundary),constraint_temporary);
	if(face->boundary) for(i = 0; i < condition_n_variables(boundary_condition(face->boundary)); i ++) n_constraints[constraint_temporary[i]]++;

	if(face->node)
	{
		printf("    face->node\n       ");
		for(i = 0; i < 2; i ++) printf(" %i",node_index(face->node[i]));
		printf("\n");
	}
	if(face->border)
	{
		printf("    face->border\n       ");
		for(i = 0; i < face->n_borders; i ++) printf(" %i",element_index(face->border[i]));
		printf("\n");
	}
	if(face->normal)
	{
		printf("    face->normal\n       ");
		for(i = 0; i < 2; i ++) printf(" %g",X_GT_EPS(face->normal[i]));
		printf("\n");
	}
	if(face->centre)
	{
		printf("    face->centre\n       ");
		for(i = 0; i < 2; i ++) printf(" %g",face->centre[i]);
		printf("\n");
	}
	if(face->size != 0.0)
		printf("    face->size\n        %g\n",face->size);
	if(face->X && face->W)
	{
		printf("    face->X");
		for(i = 0; i < face->n_quadrature; i ++) { printf("\n       "); for(j = 0; j < 2; j ++) { printf(" %+e",face->X[j][i]); } }
		printf("\n");
		printf("    face->W");
		for(i = 0; i < face->n_quadrature; i ++) printf("\n        %+e",face->W[i]);
		printf("\n");
	}
	if(face->system_index != -1)
		printf("    face->system_index\n        %i\n",face->system_index);
	if(face->Q)
	{
		printf("    face->Q");
		for(i = 0; i < solver_n_interpolations(); i ++) {
			for(j = 0; j < face->n_borders*n_bases[interpolation_variable[i]] + n_constraints[interpolation_variable[i]]*face->n_quadrature; j ++) {
				printf("\n       ");
				for(k = 0; k < face->n_quadrature; k ++) {
					printf(" %+e",X_GT_EPS(face->Q[i][j][k]));
				}
			}
			printf("\n");
		}
	}

	free(n_bases);
	free(interpolation_variable);
	free(n_constraints);
	free(constraint_temporary);
}

////////////////////////////////////////////////////////////////////////////////

void face_plot(FACE face)
{
	int i, j, k, l;
	NODE n[2];
	FACE f[ELEMENT_MAX_N_FACES];

	printf("set size ratio -1;\n");
	printf("plot '-' w lp title 'line', '-' w p title 'quadrature', '-' w p title 'centre', '-' w vec title '10%% normal'");
	for(i = 0; i < face->n_borders; i ++) printf(", '-' w l title 'border %i'",i);
	printf("\n");

	double x[2][2];
	for(i = 0; i < 2; i ++) node_x(face->node[i],x[i]);

	for(i = 0; i < 2; i ++)
	{
		for(j = 0; j < 2; j ++) printf("%e ",x[i][j]);
		printf("\n");
	}
	printf("e\n");

	for(i = 0; i < face->n_quadrature; i ++)
	{
		for(j = 0; j < 2; j ++)
		{
			printf("%e ",face->X[j][i]);
		}
		printf("\n");
	}
	printf("e\n");

	for(i = 0; i < 2; i ++) printf("%e ",face->centre[i]);
	printf("\ne\n");

	for(i = 0; i < 2; i ++) printf("%e ",face->centre[i]);
	for(i = 0; i < 2; i ++) printf("%e ",0.1*face->size*face->normal[i]);
	printf("\ne\n");

	for(i = 0; i < face->n_borders; i ++)
	{
		element_face(face->border[i],f);
		for(j = 0; j < element_n_faces(face->border[i]); j ++)
		{
			if(f[j] != face)
			{
				face_node(f[j],n);

				for(k = 0; k < 2; k ++) node_x(n[k],x[k]);

				for(k = 0; k < 2; k ++)
				{
					for(l = 0; l < 2; l ++)
					{
						printf("%e ",x[k][l]);
					}               
					printf("\n");
				}       
				printf("\n");
			}
		}               
		printf("e\n");
	}
}

////////////////////////////////////////////////////////////////////////////////

void face_free(FACE face)
{
	free(face->node);
	free(face->border);
	free(face->normal);
	free(face->centre);
	matrix_free((void *)face->X);
	free(face->W);
	tensor_free((void *)face->Q);
	free(face);
}

////////////////////////////////////////////////////////////////////////////////
