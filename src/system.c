//////////////////////////////////////////////////////////////////

#include "isitek.h"

#include "constants.h"
#include "expression.h"
#include "linear.h"

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

	/**n_u = 0;
	for(e = 0; e < n_elements; e ++)
	{
		element[e].unknown = allocate_element_unknown(&element[e],n_variables,n_basis);
		exit_if_false(element[e].unknown != NULL,"allocating element unknown indices");
	}

	for(v = 0; v < n_variables; v ++)
	{
		for(e = 0; e < n_elements; e ++)
		{
			for(i = 0; i < n_basis[v]; i ++)
			{
				if(n_u_old) if(v < n_variables_old && i < n_basis_old[v]) (*u)[*n_u+i] = u_old[element[e].unknown[v][i]];
				element[e].unknown[v][i] = *n_u+i;
			}
			*n_u += n_basis[v];
		}
	}*/

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
	int e, v, i, z = powers_taylor[0][0];

	for(e = 0; e < n_elements; e ++)
		for(v = 0; v < n_variables; v ++)
			for(i = 0; i < ORDER_TO_N_BASIS(variable_order[v]); i ++)
				u[element[e].unknown[v][i]] = initial[v]*(i == z);
}

//////////////////////////////////////////////////////////////////

void initialise_system(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, int n_u, SPARSE *system)
{
	int e, i, j, k, v;

	// allocate the system
	exit_if_false(*system = sparse_allocate(*system, n_u),"allocating the system");

	// allocate the rows
	int *row_n_non_zeros = (int *)malloc(n_u * sizeof(int)), n_non_zeros;
	exit_if_false(row_n_non_zeros != NULL,"allocating row numbers of non-zeros");

	int sum_n_basis = 0;
	for(v = 0; v < n_variables; v ++) sum_n_basis += ORDER_TO_N_BASIS(variable_order[v]);

	for(e = 0; e < n_elements; e ++)
	{
		n_non_zeros = sum_n_basis;
		for(i = 0; i < element[e].n_faces; i ++)
			n_non_zeros += sum_n_basis * (element[e].face[i]->n_borders == 2);

		for(v = 0; v < n_variables; v ++)
			for(i = 0; i < ORDER_TO_N_BASIS(variable_order[v]); i ++)
				row_n_non_zeros[element[e].unknown[v][i]] = n_non_zeros;
	}
	
	exit_if_false(*system = sparse_allocate_rows(*system, row_n_non_zeros),"allocating the rows");

	free(row_n_non_zeros);

	// set the row indices
	int *row_index = (int *)malloc(sum_n_basis * MAX_ELEMENT_N_FACES * sizeof(int));
	exit_if_false(row_index != NULL,"allocating row indices");
	for(e = 0; e < n_elements; e ++)
	{
		n_non_zeros = 0;

		for(v = 0; v < n_variables; v ++)
			for(i = 0; i < ORDER_TO_N_BASIS(variable_order[v]); i ++)
				row_index[n_non_zeros++] = element[e].unknown[v][i];

		for(i = 0; i < element[e].n_faces; i ++)
			for(j = 0; j < element[e].face[i]->n_borders; j ++)
				if(element[e].face[i]->border[j] != &element[e])
					for(v = 0; v < n_variables; v ++)
						for(k = 0; k < ORDER_TO_N_BASIS(variable_order[v]); k ++)
							row_index[n_non_zeros++] = element[e].face[i]->border[j]->unknown[v][k];

		for(v = 0; v < n_variables; v ++)
			for(i = 0; i < ORDER_TO_N_BASIS(variable_order[v]); i ++)
				sparse_set_row_indices(*system, element[e].unknown[v][i], row_index);
	}

	free(row_index);
}

//////////////////////////////////////////////////////////////////

void calculate_system(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, int n_terms, struct TERM *term, double *u_old, double *u, SPARSE system, double *residual)
{
	int a, b, d, e, i, j, k, n, p, q, t, v, x;

	// orders and numbers of bases
	int max_variable_order = 0;
	for(v = 0; v < n_variables; v ++) max_variable_order = MAX(max_variable_order,variable_order[v]);

	int *n_basis, max_n_basis = ORDER_TO_N_BASIS(max_variable_order), sum_n_basis = 0;
	exit_if_false(n_basis = (int *)malloc(n_variables * sizeof(int)),"allocating n_basis");
	for(v = 0; v < n_variables; v ++) sum_n_basis += n_basis[v] = ORDER_TO_N_BASIS(variable_order[v]);

	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order), n_points;

	// max term variables
	int max_term_n_variables = 0;
	for(t = 0; t < n_terms; t ++) max_term_n_variables = MAX(max_term_n_variables,term[t].n_variables);

	// local dense system
	double **local_value, *local_residual;
	int local_n, ldl = (1+MAX_ELEMENT_N_FACES)*sum_n_basis;
	exit_if_false(local_value = allocate_double_matrix(NULL,sum_n_basis,(1+MAX_ELEMENT_N_FACES)*sum_n_basis),"allocating local values");
	exit_if_false(local_residual = (double *)malloc(sum_n_basis * sizeof(double)),"allocating local residuals");

	// indices into the local dense system
	int *local_element, **local_adjacent;
	exit_if_false(local_element = (int *)malloc(n_variables * sizeof(int)),"allocating element local indices");
	exit_if_false(local_adjacent = allocate_integer_matrix(NULL,MAX_ELEMENT_N_FACES,n_variables),"allocating adjcent local indices");

	// working arrays
	double *basis_value, *basis_value_old, **point_value, **point_value_old, *point_term, *point_term_old;
	exit_if_false(basis_value = (double *)malloc((2*max_n_basis + MAX_FACE_N_BOUNDARIES) * sizeof(double)),"allocating basis values");
	exit_if_false(basis_value_old = (double *)malloc((2*max_n_basis + MAX_FACE_N_BOUNDARIES) * sizeof(double)),"allocating old basis values");
	exit_if_false(point_value = allocate_double_matrix(NULL,max_term_n_variables,MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer)),"allocating point values");
	exit_if_false(point_value_old = allocate_double_matrix(NULL,max_term_n_variables,MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer)),"allocating old point values");
	exit_if_false(point_term = (double *)malloc(MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer) * sizeof(double)),"allocating point term values");
	exit_if_false(point_term_old = (double *)malloc(MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer) * sizeof(double)),"allocating old point term values");
	int lda = MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer);
	double **A = allocate_double_matrix(NULL,max_n_basis,lda);

	// expression evaluation working memory
	int expression_max_recusions = 0;
	for(t = 0; t < n_terms; t ++) expression_max_recusions = MAX(expression_max_recusions,expression_number_of_recursions(term[t].residual));
	double **expression_work;
	exit_if_false(expression_work = allocate_double_matrix(NULL,expression_max_recusions,MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer)),"allocating expression work");

	// blas parameters
	char trans[2] = "NT";
	int unit = 1;
	double one = 1.0, zero = 0.0;

	// loop over the elements
	for(e = 0; e < n_elements; e ++)
	{
		// local system indices
		local_n = 0;
		for(v = 0; v < n_variables; v ++)
		{
			local_element[v] = local_n;
			local_n += n_basis[v];
		}
		for(a = 0; a < element[e].n_faces; a ++)
			for(i = 0; i < element[e].face[a]->n_borders; i ++)
				if(element[e].face[a]->border[i] != &element[e])
					for(v = 0; v < n_variables; v ++)
					{
						local_adjacent[a][v] = local_n;
						local_n += n_basis[v];
					}

		// initialise the local system
		for(i = 0; i < sum_n_basis; i ++)
			for(j = 0; j < local_n; j ++)
				local_value[i][j] = 0.0;
		for(i = 0; i < sum_n_basis; i ++)
			local_residual[i] = 0.0;

		// integration points
		n_points = (element[e].n_faces - 2)*n_hammer;

		// loop over the terms
		for(t = 0; t < n_terms; t ++)
		{
			// equation and direction
			q = term[t].equation;
			x = powers_taylor[term[t].type == 'x'][term[t].type == 'y'];

			// FEM component -----------------------//

			// values at the intergration points
			for(i = 0; i < term[t].n_variables; i ++)
			{
				v = term[t].variable[i]; d = term[t].differential[i];

				for(j = 0; j < n_basis[v]; j ++)
				{
					basis_value[j] = u[element[e].unknown[v][j]];
					basis_value_old[j] = u_old[element[e].unknown[v][j]];
				}

				dgemv_(&trans[0],&n_points,&n_basis[v],
						&one,
						element[e].P[d][0],&n_points,
						basis_value,&unit,
						&zero,
						point_value[i],&unit);
				dgemv_(&trans[0],&n_points,&n_basis[v],
						&one,
						element[e].P[d][0],&n_points,
						basis_value_old,&unit,
						&zero,
						point_value_old[i],&unit);
			}

			// jacobian
			for(i = 0; i < term[t].n_variables; i ++)
			{
				v = term[t].variable[i]; d = term[t].differential[i];

				expression_evaluate(n_points, point_term, term[t].jacobian[i], point_value, expression_work);
				for(p = 0; p < n_points; p ++) point_term[p] *= - element[e].W[p] * term[t].implicit;

				for(j = 0; j < n_basis[q]; j ++)
					for(p = 0; p < n_points; p ++)
						A[j][p] = element[e].P[x][j][p] * point_term[p];

				dgemm_(&trans[1],&trans[0],&n_basis[v],&n_basis[q],&n_points,
						&one,
						element[e].P[d][0],&n_points,
						A[0],&lda,
						&one,
						&local_value[local_element[q]][local_element[v]],&ldl);
			}

			// residual
			expression_evaluate(n_points, point_term, term[t].residual, point_value, expression_work);
			expression_evaluate(n_points, point_term_old, term[t].residual, point_value_old, expression_work);
			for(p = 0; p < n_points; p ++) point_term[p] = - element[e].W[p] *
				(term[t].implicit * point_term[p] + (1.0 - term[t].implicit) * point_term_old[p]);
			dgemv_(&trans[1],&n_points,&n_basis[q],
					&one,
					element[e].P[x][0],&n_points,
					point_term,&unit,
					&one,
					&local_residual[local_element[q]],&unit);
			
			// break if a source term
			if(term[t].type == 's') continue;

			// DG FEM component --------------------//

			// direction
			x = term[t].type == 'y';

			// loop over the adjacent elements
			for(a = 0; a < element[e].n_faces; a ++)
			{
				// values at the intergration points
				for(i = 0; i < term[t].n_variables; i ++)
				{
					v = term[t].variable[i]; d = term[t].differential[i];

					for(j = 0; j < element[e].face[a]->n_borders; j ++)
					{
						for(k = 0; k < n_basis[v]; k ++)
						{
							basis_value[k+j*n_basis[v]] = u[element[e].face[a]->border[j]->unknown[v][k]];
							basis_value_old[k+j*n_basis[v]] = u_old[element[e].face[a]->border[j]->unknown[v][k]];
						}
					}
					for(j = 0; j < element[e].face[a]->n_boundaries[v]; j ++)
					{
						basis_value[j+element[e].face[a]->n_borders*n_basis[v]] = element[e].face[a]->boundary[v][j]->value;
						basis_value_old[j+element[e].face[a]->n_borders*n_basis[v]] = element[e].face[a]->boundary[v][j]->value;
					}

					if(term[t].method[i] == 'r' || d != powers_taylor[0][0] || element[e].face[a]->n_borders == 0 || element[e].face[a]->n_boundaries[v])
					{
						n = element[e].face[a]->n_borders*n_basis[v] + element[e].face[a]->n_boundaries[v];
						dgemv_(&trans[0],&n_gauss,&n,
								&one,
								element[e].face[a]->Q[v][d][0],&n_gauss,
								basis_value,&unit,
								&zero,
								point_value[i],&unit);
						dgemv_(&trans[0],&n_gauss,&n,
								&one,
								element[e].face[a]->Q[v][d][0],&n_gauss,
								basis_value_old,&unit,
								&zero,
								point_value_old[i],&unit);
					}
				}

				// jacobian
				for(i = 0; i < term[t].n_variables; i ++)
				{
					v = term[t].variable[i]; d = term[t].differential[i];

					expression_evaluate(n_gauss, point_term, term[t].jacobian[i], point_value, expression_work);
					for(p = 0; p < n_gauss; p ++) point_term[p] *= element[e].orient[a] *
						element[e].face[a]->normal[x] * term[t].implicit * element[e].face[a]->W[p];

					for(j = 0; j < n_basis[q]; j ++)
						for(p = 0; p < n_gauss; p ++)
							A[j][p] = element[e].Q[a][j][p] * point_term[p];

					if(term[t].method[i] == 'r' || d != powers_taylor[0][0] || element[e].face[a]->n_borders == 0 || element[e].face[a]->n_boundaries[v])
					{
						for(j = 0; j < element[e].face[a]->n_borders; j ++)
						{
							b = element[e].face[a]->border[j] == &element[e] ? local_element[v] : local_adjacent[a][v];
							dgemm_(&trans[1],&trans[0],&n_basis[v],&n_basis[q],&n_gauss,
									&one,
									element[e].face[a]->Q[v][d][j*n_basis[v]],&n_gauss,
									A[0],&lda,
									&one,
									&local_value[local_element[q]][b],&ldl);
						}
					}
				}

				// residual
				expression_evaluate(n_gauss, point_term, term[t].residual, point_value, expression_work);
				expression_evaluate(n_gauss, point_term_old, term[t].residual, point_value_old, expression_work);
				for(p = 0; p < n_gauss; p ++) point_term[p] = element[e].orient[a] * element[e].face[a]->normal[x] * element[e].face[a]->W[p] *
					(term[t].implicit * point_term[p] + (1.0 - term[t].implicit) * point_term_old[p]);
				dgemv_(&trans[1],&n_gauss,&n_basis[q],
						&one,
						element[e].Q[a][0],&n_gauss,
						point_term,&unit,
						&one,
						&local_residual[local_element[q]],&unit);
			}
		}

		// add to the global system
		for(v = 0; v < n_variables; v ++)
		{
			for(i = 0; i < n_basis[v]; i ++)
			{
				sparse_set_row_values(system, element[e].unknown[v][i], local_value[local_element[v]+i]);
				residual[element[e].unknown[v][i]] = local_residual[local_element[v]+i];
			}
		}
	}

	free(n_basis);

	free(local_residual);
	destroy_matrix((void *)local_value);
	free(local_element);
	destroy_matrix((void *)local_adjacent);

	free(basis_value);
	free(basis_value_old);
	destroy_matrix((void *)point_value);
	destroy_matrix((void *)point_value_old);
	free(point_term);
	free(point_term_old);
	destroy_matrix((void *)A);

	destroy_matrix((void *)expression_work);
}

//////////////////////////////////////////////////////////////////

void calculate_maximum_residuals(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, double *residual, double *max_residual)
{
	int e, i, v;

	for(v = 0; v < n_variables; v ++) max_residual[v] = 0.0;

	for(e = 0; e < n_elements; e ++)
		for(v = 0; v < n_variables; v ++)
			for(i = 0; i < ORDER_TO_N_BASIS(variable_order[v]); i ++)
				max_residual[v] = MAX(max_residual[v],fabs(residual[element[e].unknown[v][i]]));
}

//////////////////////////////////////////////////////////////////
