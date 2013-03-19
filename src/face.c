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

const NODE * face_node(FACE face)
{
	return face->node;
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

const ELEMENT * face_border(FACE face)
{
	return face->border;
}

////////////////////////////////////////////////////////////////////////////////

void face_calculate_size(FACE face)
{
	int i;
	double dx;
	face->size = 0;
	for(i = 0; i < 2; i ++)
	{
		dx = node_x(face->node[1])[i] - node_x(face->node[0])[i];
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
	for(i = 0; i < 2; i ++)
	{
		face->centre[i] = 0.5*(node_x(face->node[1])[i] + node_x(face->node[0])[i]);
	}

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

const double * face_centre(FACE face)
{
	return face->centre;
}

////////////////////////////////////////////////////////////////////////////////

int face_calculate_normal(FACE face)
{
	face->normal = (double *)malloc(2*sizeof(double));
	if(face->normal == NULL) return FACE_MEMORY_ERROR;

	face->normal[0] = (- node_x(face->node[1])[1] + node_x(face->node[0])[1])/face->size;
	face->normal[1] = (+ node_x(face->node[1])[0] - node_x(face->node[0])[0])/face->size;

	return FACE_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

const double * face_normal(FACE face)
{
	return face->normal;
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

	for(i = 0; i < solver_n_gauss(); i ++)
	{
		for(j = 0; j < 2; j ++)
		{
			face->X[j][i] =
				0.5 * (1.0 - quadrature_gauss_location(solver_n_gauss()-1,i)) * node_x(face->node[0])[j] +
				0.5 * (1.0 + quadrature_gauss_location(solver_n_gauss()-1,i)) * node_x(face->node[1])[j];
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

double const * const * face_quadrature_x(FACE face)
{
	return face->X;
}

////////////////////////////////////////////////////////////////////////////////

const double * face_quadrature_w(FACE face)
{
	return face->W;
}

////////////////////////////////////////////////////////////////////////////////

int face_add_to_system(FACE face, SPARSE system)
{
	int i, j;

	int **unknown = matrix_integer_new(NULL,face->n_borders,solver_variable_sum_n_bases());
	if(unknown == NULL) return FACE_MEMORY_ERROR;

	for(i = 0; i < face->n_borders; i ++)
		for(j = 0; j < solver_variable_sum_n_bases(); j ++)
			unknown[i][j] = element_unknown(face->border[i])[0][j];

	face->system_index = sparse_insert_sub_matrix(system,face->n_borders*solver_variable_sum_n_bases(),unknown[0],face->n_borders*solver_variable_sum_n_bases(),unknown[0]);

	matrix_free((void *)unknown);

	return face->system_index >= 0 ? FACE_SUCCESS : FACE_LOGIC_ERROR;
}

////////////////////////////////////////////////////////////////////////////////

void face_print(FACE face)
{
	printf("face %i\n",face->index);

	int i, j, k;

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
		for(i = 0; i < solver_n_interpolations(); i ++)
		{
			int n = face->n_borders*solver_variable_n_bases()[solver_interpolation_variable()[i]];
			if(face->boundary) n += condition_variable_n_constraints(boundary_condition(face->boundary))[solver_interpolation_variable()[i]]*face->n_quadrature;

			for(j = 0; j < n; j ++)
			{
				printf("\n       ");
				for(k = 0; k < face->n_quadrature; k ++)
				{
					printf(" %+e",X_GT_EPS(face->Q[i][j][k]));
				}
			}
			printf("\n");
		}
	}
}

////////////////////////////////////////////////////////////////////////////////

void face_plot(FACE face)
{
	int i, j, k, l;

	printf("set size ratio -1;\n");
	printf("plot '-' w lp title 'line', '-' w p title 'quadrature', '-' w p title 'centre', '-' w vec title '10%% normal'");
	for(i = 0; i < face->n_borders; i ++) printf(", '-' w l title 'border %i'",i);
	printf("\n");

	for(i = 0; i < 2; i ++)
	{
		for(j = 0; j < 2; j ++) printf("%e ",node_x(face->node[i])[j]);
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
		for(j = 0; j < element_n_faces(face->border[i]); j ++)
		{
			if(element_face(face->border[i])[j] != face)
			{
				for(k = 0; k < 2; k ++)
				{
					for(l = 0; l < 2; l ++)
					{
						printf("%e ",node_x(face_node(element_face(face->border[i])[j])[k])[l]);
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
