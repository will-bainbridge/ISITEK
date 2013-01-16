////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "element.h"
#include "memory.h"
#include "numerics.h"
#include "quadrature.h"
#include "solver.h"

////////////////////////////////////////////////////////////////////////////////

#include "element.r"

////////////////////////////////////////////////////////////////////////////////

ELEMENT element_new(int index)
{
	ELEMENT element = (ELEMENT)malloc(sizeof(struct s_ELEMENT));
	if(element == NULL) return NULL;

	element->index = index;

	element->n_faces = 0;
	element->face = NULL;
	
	element->centre = NULL;
	element->size = 0.0;

	element->n_quadrature = 0;
	element->X = NULL;
	element->W = NULL;

	element->unknown = NULL;

	element->system_index = -1;

	element->P = NULL;
	element->Q = NULL;
	element->I = NULL;
	element->V = NULL;
	element->L = NULL;

	return element;
}

////////////////////////////////////////////////////////////////////////////////

int element_index(ELEMENT element)
{
	return element->index;
}

////////////////////////////////////////////////////////////////////////////////

int element_read_faces(FILE *file, FACE *face, ELEMENT element)
{
	if(fscanf(file,"%i",&(element->n_faces)) != 1)
	{
		return ELEMENT_READ_ERROR;
	}

	element->face = (FACE *)malloc(element->n_faces * sizeof(FACE));
	if(element->face == NULL) return ELEMENT_MEMORY_ERROR;

	int i, index;

	for(i = 0; i < element->n_faces; i ++)
	{
		if(fscanf(file,"%i",&index) != 1)
		{
			return ELEMENT_READ_ERROR;
		}
		element->face[i] = face[index];
	}

	return ELEMENT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int element_n_faces(ELEMENT element)
{
	return element->n_faces;
}

////////////////////////////////////////////////////////////////////////////////

void element_face(ELEMENT element, FACE * face)
{
	int i;
	for(i = 0; i < element->n_faces; i ++) face[i] = element->face[i];
}

////////////////////////////////////////////////////////////////////////////////

int element_add_border(ELEMENT element)
{
	// counters
	int i, j;

	// get an ordered list of faces and orientations around the element 
	NODE n[2], n_old[2];
	FACE f[element->n_faces];
	int o[element->n_faces];
	f[0] = element->face[0];
	o[0] = 0;
	for(i = 1; i < element->n_faces; i ++)
	{
		for(j = 0; j < element->n_faces; j ++)
		{
			f[i] = element->face[j];
			if(f[i] == f[i-1]) continue;
			face_node(f[i-1],n_old);
			face_node(f[i],n);
			if(n_old[o[i-1]] == n[0]) { o[i] = 1; break; }
			if(n_old[o[i-1]] == n[1]) { o[i] = 0; break; }
		}
	}

	// check area and flip all orientations if -ve
	double a = 0;
	for(i = 0; i < element->n_faces; i ++)
	{
		face_node(f[i],n);
		a += node_x(n[!o[i]],0)*node_x(n[o[i]],1) - node_x(n[o[i]],0)*node_x(n[!o[i]],1);
	}
	if(a < 0) for(i = 0; i < element->n_faces; i ++) o[i] = !o[i];

	// add face borders
	for(i = 0; i < element->n_faces; i ++)
	{
		if(face_add_border(f[i],element,o[i]) != FACE_SUCCESS)
		{
			return ELEMENT_FAIL;
		}
	}

	return ELEMENT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int element_calculate_quadrature(ELEMENT element)
{
	int i, j, k;
	double x[2], dx[2][2], size;
	NODE n[3], temp;
	FACE *f = &element->face[0];
	ELEMENT e[2];

	// allocation
	element->n_quadrature = (element->n_faces-2)*solver_n_hammer();
	element->X = matrix_double_new(element->X, 2, element->n_quadrature);
	if(element->X == NULL) return ELEMENT_MEMORY_ERROR;
	element->W = (double *)malloc(element->n_quadrature * sizeof(double));
	if(element->W == NULL) return ELEMENT_MEMORY_ERROR;

	// triangulation base point
	face_node(*f,n);
	for(i = 0; i < 2; i ++) x[i] = node_x(n[0],i);

	// opposite edges
	i = 0;
	while(i < element->n_faces - 2)
	{
		f += 1;

		// continue if any nodes match the base point
		face_node(*f,&n[1]);
		face_border(*f,e);
		if(n[1] == n[0] || n[2] == n[0]) continue;
		if(e[0] == element) { temp = n[1]; n[1] = n[2]; n[2] = temp; }

		// edge vectors
		for(j = 0; j < 2; j ++)
			for(k = 0; k < 2; k ++)
				dx[j][k] = node_x(n[j+1],k) - x[k];

		// size
		size = dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0];

		// hammer points and weights
		for(j = 0; j < solver_n_hammer(); j ++)
		{
			for(k = 0; k < 2; k ++)
			{
				element->X[k][i*solver_n_hammer()+j] = x[k] +
					quadrature_hammer_location(solver_n_hammer()-1,0,j)*dx[0][k] +
					quadrature_hammer_location(solver_n_hammer()-1,1,j)*dx[1][k];
			}
			element->W[i*solver_n_hammer()+j] = quadrature_hammer_weight(solver_n_hammer()-1,j) * size;
		}

		// next edge
		i ++;
	}

	return ELEMENT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int element_n_quadrature(ELEMENT element)
{
        return element->n_quadrature;
}

////////////////////////////////////////////////////////////////////////////////

void element_quadrature_x(ELEMENT element, double **x)
{
        int i, j;
        for(i = 0; i < element->n_quadrature; i ++)
                for(j = 0; j < 2; j ++)
                        x[j][i] = element->X[j][i];
}

////////////////////////////////////////////////////////////////////////////////

void element_quadrature_w(ELEMENT element, double *w)
{
        int i;
        for(i = 0; i < element->n_quadrature; i ++) w[i] = element->W[i];
}

////////////////////////////////////////////////////////////////////////////////

void element_calculate_size(ELEMENT element)
{
	int i;
	element->size = 0;
	for(i = 0; i < element->n_quadrature; i ++) element->size += element->W[i];
	element->size = sqrt(element->size);
}

////////////////////////////////////////////////////////////////////////////////

double element_size(ELEMENT element)
{
	return element->size;
}

////////////////////////////////////////////////////////////////////////////////

int element_calculate_centre(ELEMENT element)
{
	element->centre = (double *)malloc(2*sizeof(double));
	if(element->centre == NULL) return ELEMENT_MEMORY_ERROR;

	int i, j;
	for(i = 0; i < 2; i ++) element->centre[i] = 0;
	for(i = 0; i < element->n_quadrature; i ++) for(j = 0; j < 2; j ++) element->centre[j] += element->W[i] * element->X[j][i];
	for(i = 0; i < 2; i ++) element->centre[i] /= element->size*element->size;

	return ELEMENT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

void element_centre(ELEMENT element, double *x)
{
	int i;
	for(i = 0; i < 2; i ++) x[i] = element->centre[i];
}

////////////////////////////////////////////////////////////////////////////////

int element_calculate_unknowns(ELEMENT element, int n_elements)
{
	int i, j, offset = 0;

	int n_variables = solver_n_variables(), *n_bases = (int *)malloc(n_variables * sizeof(int)), sum_n_bases = solver_variable_sum_n_bases();
	if(n_bases == NULL) return ELEMENT_MEMORY_ERROR;
	solver_variable_n_bases(n_bases);

	element->unknown = (int **)malloc(n_variables * sizeof(int *));
	if(element->unknown == NULL) return ELEMENT_MEMORY_ERROR;
	element->unknown[0] = (int *)malloc(sum_n_bases * sizeof(int));
	if(element->unknown[0] == NULL) return ELEMENT_MEMORY_ERROR;
	for(i = 1; i < n_variables; i++) element->unknown[i] = element->unknown[i-1] + n_bases[i-1];

	// number by element then variable
	for(i = 0; i < n_variables; i ++)
	{
		for(j = 0; j < n_bases[i]; j ++)
		{
			element->unknown[i][j] = offset + element->index*n_bases[i] + j;
		}
		offset += n_elements*n_bases[i];
	}

	// number by variable then element
	/*for(i = 0; i < n_variables; i ++)
	{
		for(j = 0; j < n_bases[i]; j ++)
		{
			element->unknown[i][j] = element->index*sum_n_bases + offset + j;
		}
		offset += n_bases[i];
	}*/

	free(n_bases);

	return ELEMENT_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

void element_unknown(ELEMENT element, int **unknown)
{
	int i, j, m = solver_n_variables(), n;

	solver_variable_n_bases(unknown[0]);

	for(i = m - 1; i > 0; i --) unknown[i][0] = unknown[0][i];

	for(i = 0; i < m; i ++)
	{
		n = unknown[i][0];
		for(j = 0; j < n; j ++) unknown[i][j] = element->unknown[i][j];
	}
}

////////////////////////////////////////////////////////////////////////////////

int element_add_to_system(ELEMENT element, SPARSE system)
{
	element->system_index = sparse_insert_sub_matrix(system,solver_variable_sum_n_bases(),element->unknown[0],solver_variable_sum_n_bases(),element->unknown[0]);

	return element->system_index >= 0 ? ELEMENT_SUCCESS : ELEMENT_LOGIC_ERROR;
}

////////////////////////////////////////////////////////////////////////////////

void element_print(ELEMENT element)
{
	printf("element %i\n",element->index);

	int *n_bases = (int *)malloc(solver_n_variables() * sizeof(int));
	solver_variable_n_bases(n_bases);

	int i, j, k;
	if(element->face)
	{
		printf("    element->face\n       ");
		for(i = 0; i < element->n_faces; i ++) printf(" %i",face_index(element->face[i]));
		printf("\n");
	}
	if(element->X && element->W)
	{
		printf("    element->X");
		for(i = 0; i < element->n_quadrature; i ++) { printf("\n       "); for(j = 0; j < 2; j ++) printf(" %+e",element->X[j][i]); }
		printf("\n");
		printf("    element->W");
		for(i = 0; i < element->n_quadrature; i ++) printf("\n        %+e",element->W[i]);
		printf("\n");
	}
	if(element->size != 0.0)
		printf("    element->size\n        %g\n",element->size);
	if(element->centre)
	{
		printf("    element->centre\n       ");
		for(i = 0; i < 2; i ++) printf(" %g",element->centre[i]);
		printf("\n");
	}
	if(element->unknown)
	{
		printf("    element->unknown");
		for(i = 0; i < solver_n_variables(); i ++) { printf("\n       "); for(j = 0; j < n_bases[i]; j ++) printf(" %i",element->unknown[i][j]); }
		printf("\n");
	}
	if(element->system_index != -1)
		printf("    element->system_index\n        %i\n",element->system_index);
	if(element->P)
	{
		printf("    element->P");
		for(i = 0; i < solver_variable_max_n_bases(); i ++) {
			for(j = 0; j < solver_variable_max_n_bases(); j ++) {
				printf("\n       ");
				for(k = 0; k < element->n_quadrature; k ++) {
					printf(" %+e",element->P[i][j][k]);
				}
			}
			printf("\n");
		}
	}
	if(element->Q)
	{
		printf("    element->Q");
		for(i = 0; i < element->n_faces; i ++) {
			for(j = 0; j < solver_variable_max_n_bases(); j ++) {
				printf("\n       ");
				for(k = 0; k < face_n_quadrature(element->face[i]); k ++) {
					printf(" %+e",element->Q[i][j][k]);
				}
			}
			printf("\n");
		}
	}
	if(element->I)
	{
		printf("    element->I");
		for(i = 0; i < solver_n_variables(); i ++) {
			for(j = 0; j < element->n_quadrature; j ++) {
				printf("\n       ");
				for(k = 0; k < n_bases[i]; k ++) {
					printf(" %+e",element->I[i][j][k]);
				}
			}
			printf("\n");
		}
	}
	if(element->V)
	{
		printf("    element->V");
		for(i = 0; i < solver_variable_max_n_bases(); i ++) {
			printf("\n       ");
			for(j = 0; j < element->n_faces; j ++) {
				printf(" %+e",element->V[i][j]);
			}
		}
		printf("\n");
	}

	if(element->L)
	{
		printf("    element->L");
		for(i = 0; i < solver_n_variables(); i ++) {
			for(j = 0; j < n_bases[i]; j ++) {
				printf("\n       ");
				for(k = 0; k < n_bases[i]; k ++) {
					printf(" %+e",element->L[i][j][k]);
				}
			}
			printf("\n");
		}
	}

	free(n_bases);
}

////////////////////////////////////////////////////////////////////////////////

void element_plot(ELEMENT element)
{
	int i, j, k;
	NODE node[2];

	printf("set size ratio -1;\n");
	printf("plot '-' w lp title 'boundary', '-' w p title 'quadrature', '-' w p title 'centre'\n");

	for(i = 0; i < element->n_faces; i ++)
	{
		face_node(element->face[i],node);

		for(j = 0; j < 2; j ++)
		{
			for(k = 0; k < 2; k ++)
			{
				printf("%e ",node_x(node[j],k));
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("e\n");

	for(i = 0; i < element->n_quadrature; i ++)
	{
		for(j = 0; j < 2; j ++)
		{
			printf("%e ",element->X[j][i]);
		}
		printf("\n");
	}
	printf("e\n");

	for(i = 0; i < 2; i ++) printf("%e ",element->centre[i]);
	printf("\ne\n");
}

////////////////////////////////////////////////////////////////////////////////

void element_free(ELEMENT element)
{
	free(element->face);
	free(element->centre);
	matrix_free((void *)element->X);
	free(element->W);
	if(element->unknown) free(element->unknown[0]);
	free(element->unknown);
	tensor_free((void *)element->P);
	tensor_free((void *)element->Q);
	matrix_free((void *)element->V);
	tensor_free((void *)element->I);
	tensor_free((void *)element->L);
	free(element);
}

////////////////////////////////////////////////////////////////////////////////
