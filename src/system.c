//////////////////////////////////////////////////////////////////

#include "isitek.h"

//////////////////////////////////////////////////////////////////

void update_element_unknowns(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element, int n_u_old, int *n_u, double *u_old, double **u)
{
	int e,v,i;

	int *n_basis = (int *)malloc(n_variables * sizeof(int));
	exit_if_false(n_basis != NULL,"allocating n_basis");
	for(v = 0; v < n_variables; v ++) n_basis[v] = ORDER_TO_N_BASIS(variable_order[v]);

	*n_u = 0;
	for(e = 0; e < n_elements; e ++) for(v = 0; v < n_variables; v ++) *n_u += n_basis[v];
	*u = (double *)malloc(*n_u * sizeof(double));
	exit_if_false(*u != NULL,"allocating u");

	int *n_basis_old = (int *)malloc(n_variables * sizeof(int));
	exit_if_false(n_basis_old != NULL,"allocating n_basis_old");
	for(v = 0; v < n_variables; v ++) n_basis_old[v] = 0;
	for(v = 0; v < MIN(n_variables,n_variables_old); v ++) n_basis_old[v] = ORDER_TO_N_BASIS(variable_order_old[v]);

	*n_u = 0;
	for(e = 0; e < n_elements; e ++)
	{
		element[e].unknown = allocate_element_unknown(&element[e],n_variables,n_basis);
		exit_if_false(element[e].unknown != NULL,"allocating element unknown indices");

		for(v = 0; v < n_variables; v ++)
		{
			for(i = 0; i < n_basis[v]; i ++)
			{
				if(n_u_old) if(v < n_variables_old && i < n_basis_old[v]) (*u)[*n_u+i] = u_old[element[e].unknown[v][i]];
				element[e].unknown[v][i] = *n_u+i;
			}
			*n_u += n_basis[v];
		}
	}

	/*e = 10;

	for(v = 0; v < n_variables; v ++){
		for(i = 0; i < n_basis[v]; i ++){
			printf("%i ",element[e].unknown[v][i]);
		} printf("\n");
	} printf("\n");

	variable_order[1] = 3;
	n_basis[1] = ORDER_TO_N_BASIS(variable_order[1]);

	element[e].unknown = allocate_element_unknown(&element[e],n_variables,n_basis);

	for(v = 0; v < n_variables; v ++){
		for(i = 0; i < n_basis[v]; i ++){
			printf("%i ",element[e].unknown[v][i]);
		} printf("\n");
	} printf("\n");*/

	free(n_basis);
	free(n_basis_old);
}

//////////////////////////////////////////////////////////////////

void update_face_boundaries(int n_variables, int n_faces, struct FACE *face, int n_boundaries, struct BOUNDARY *boundary)
{
	int b, f, i, v;

	for(f = 0; f < n_faces; f ++)
	{
		exit_if_false(face[f].n_boundaries = allocate_face_n_boundaries(&face[f],n_variables),"allocating face numbers of boundaries");
		for(v = 0; v < n_variables; v ++) face[f].n_boundaries[v] = 0;
		exit_if_false(face[f].boundary = allocate_face_boundary(&face[f],n_variables),"allocating face boundaries");
	}

	struct FACE *p;
	
	for(b = 0; b < n_boundaries; b ++)
	{
		v = boundary[b].variable;

		for(i = 0; i < boundary[b].n_faces; i ++)
		{
			p = boundary[b].face[i];

			p->n_boundaries[v] ++;

			exit_if_false(p->boundary = allocate_face_boundary(p,n_variables),"allocating face boundaries");

			p->boundary[v][p->n_boundaries[v] - 1] = &boundary[b];
		}
	}

	/*int n;
	for(f = 0; f < n_faces; f ++)
	{
		n = 0;
		for(v = 0; v < n_variables; v ++) n += face[f].n_boundaries[v];

		if(n)
		{
			printf("%i ",f);
			for(v = 0; v < n_variables; v ++)
			{
				printf("( ");
				for(i = 0; i < face[f].n_boundaries[v]; i ++) printf("%li ",face[f].boundary[v][i] - &boundary[0]);
				printf(") ");
			}
			printf("\n");
		}
	}*/
}

//////////////////////////////////////////////////////////////////

void initialise_values(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, double *initial, double *u)
{
	int e, v, i;

	for(e = 0; e < n_elements; e ++)
		for(v = 0; v < n_variables; v ++)
			for(i = 0; i < ORDER_TO_N_BASIS(variable_order[v]); i ++)
				u[element[e].unknown[v][i]] = initial[v];
}

//////////////////////////////////////////////////////////////////
