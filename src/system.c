//////////////////////////////////////////////////////////////////

#include "isitek.h"

#include "constants.h"
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

void initialise_values(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, EXPRESSION *initial, double *u)
{
	int e, i, v;

	int max_variable_order = 0;
	for(v = 0; v < n_variables; v ++) max_variable_order = MAX(max_variable_order,variable_order[v]);
	int *n_basis, max_n_basis = ORDER_TO_N_BASIS(max_variable_order), sum_n_basis = 0;
	exit_if_false(n_basis = (int *)malloc(n_variables * sizeof(int)),"allocating n_basis");
	for(v = 0; v < n_variables; v ++) sum_n_basis += n_basis[v] = ORDER_TO_N_BASIS(variable_order[v]);
	int n_hammer = ORDER_TO_N_HAMMER(max_variable_order), n_points;

	int expression_max_recusions = 1;
	for(v = 0; v < n_variables; v ++) expression_max_recusions = MAX(expression_max_recusions,expression_number_of_recursions(initial[v]));
	double **expression_work;
	exit_if_false(expression_work = allocate_double_matrix(NULL,expression_max_recusions,(MAX_ELEMENT_N_FACES-2)*n_hammer),"allocating expression work");

	double *point_initial = (double *)malloc((MAX_ELEMENT_N_FACES-2)*n_hammer * sizeof(double));
	double *basis_initial = (double *)malloc(max_n_basis * sizeof(double));

	char trans[2] = "NT";
	int int_1 = 1;
	double dbl_0 = 0.0, dbl_1 = 1.0;

	for(e = 0; e < n_elements; e ++)
	{
		n_points = n_hammer * (element[e].n_faces - 2);

		for(v = 0; v < n_variables; v ++)
		{
			expression_evaluate(n_points,point_initial,initial[v],element[e].X,expression_work);
			dgemv_(&trans[0],&n_basis[v],&n_points,
					&dbl_1,
					element[e].I[v][0],&n_basis[v],
					point_initial,&int_1,
					&dbl_0,
					basis_initial,&int_1);

			for(i = 0; i < n_basis[v]; i ++) u[element[e].unknown[v][i]] = basis_initial[i];
		}
	}

	free(n_basis);
	destroy_matrix((void *)expression_work);
	free(point_initial);
	free(basis_initial);
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

void calculate_system(int n_variables, int *variable_order, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element, int n_terms, struct TERM *term, int n_u, double *u_old, double *u, SPARSE system, double *residual)
{
	int b, d, e, f, i, j, k, n, p, q, t, v, x;

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
	int max_point_n_variables = max_term_n_variables + EXPRESSION_VARIABLE_INDEX;

	// local dense system
	double ***local_value, **local_residual;
	int local_n, ldl = (1+MAX_ELEMENT_N_FACES)*sum_n_basis;
	exit_if_false(local_value = allocate_double_tensor(NULL,MAX_FACE_N_BORDERS,sum_n_basis,(1+MAX_ELEMENT_N_FACES)*sum_n_basis),"allocating local values");
	exit_if_false(local_residual = allocate_double_matrix(NULL,MAX_FACE_N_BORDERS,sum_n_basis),"allocating local residuals");

	// indices into the local dense system
	int *local_element, **local_adjacent;
	exit_if_false(local_element = (int *)malloc(n_variables * sizeof(int)),"allocating element local indices");
	exit_if_false(local_adjacent = allocate_integer_matrix(NULL,MAX_FACE_N_BORDERS,n_variables),"allocating adjacent local indices");

	// working arrays
	double *basis_value, *basis_value_old, **point_value, **point_value_old, *point_term, *point_term_old;
	exit_if_false(basis_value = (double *)malloc((2*max_n_basis + MAX_FACE_N_BOUNDARIES) * sizeof(double)),"allocating basis values");
	exit_if_false(basis_value_old = (double *)malloc((2*max_n_basis + MAX_FACE_N_BOUNDARIES) * sizeof(double)),"allocating old basis values");
	exit_if_false(point_value = allocate_double_matrix(NULL,max_point_n_variables,MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer)),"allocating point values");
	exit_if_false(point_value_old = allocate_double_matrix(NULL,max_point_n_variables,MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer)),"allocating old point values");
	exit_if_false(point_term = (double *)malloc(MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer) * sizeof(double)),"allocating point term values");
	exit_if_false(point_term_old = (double *)malloc(MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer) * sizeof(double)),"allocating old point term values");
	int lda = MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer);
	double **A = allocate_double_matrix(NULL,max_n_basis,lda);

	// expression evaluation working memory
	int expression_max_recusions = 1;
	for(t = 0; t < n_terms; t ++) expression_max_recusions = MAX(expression_max_recusions,expression_number_of_recursions(term[t].residual));
	double **expression_work;
	exit_if_false(expression_work = allocate_double_matrix(NULL,expression_max_recusions,MAX(n_gauss,(MAX_ELEMENT_N_FACES-1)*n_hammer)),"allocating expression work");

	// opposite face index
	int opposite[MAX_FACE_N_BORDERS];

	// blas parameters
	char trans[2] = "NT";
	int int_0 = 0, int_1 = 1;
	double dbl_0 = 0.0, dbl_1 = 1.0, alpha;

	// zero the system
	sparse_set_zero(system);
	dcopy_(&n_u,&dbl_0,&int_0,residual,&int_1);

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

		// zero the local system
		n = sum_n_basis;
		dcopy_(&n,&dbl_0,&int_0,local_residual[0],&int_1);
		n = sum_n_basis*ldl;
		dcopy_(&n,&dbl_0,&int_0,local_value[0][0],&int_1);

		// integration points
		n_points = (element[e].n_faces - 2)*n_hammer;

		// location values at the integration points
		for(i = 0; i < 2; i ++)
		{
			dcopy_(&n_points,element[e].X[i],&int_1,point_value[i+EXPRESSION_LOCATION_INDEX],&int_1);
			dcopy_(&n_points,element[e].X[i],&int_1,point_value_old[i+EXPRESSION_LOCATION_INDEX],&int_1);
		}

		// loop over the terms
		for(t = 0; t < n_terms; t ++)
		{
			// equation and direction
			q = term[t].equation;
			x = powers_taylor[term[t].type == 'x'][term[t].type == 'y'];

			// values at the intergration points
			for(i = 0; i < term[t].n_variables; i ++)
			{
				v = term[t].variable[i];
				d = term[t].differential[i];

				for(j = 0; j < n_basis[v]; j ++)
				{
					basis_value[j] = u[element[e].unknown[v][j]];
					basis_value_old[j] = u_old[element[e].unknown[v][j]];
				}

				dgemv_(&trans[0],&n_points,&n_basis[v],
						&dbl_1,
						element[e].P[d][0],&n_points,
						basis_value,&int_1,
						&dbl_0,
						point_value[i+EXPRESSION_VARIABLE_INDEX],&int_1);
				dgemv_(&trans[0],&n_points,&n_basis[v],
						&dbl_1,
						element[e].P[d][0],&n_points,
						basis_value_old,&int_1,
						&dbl_0,
						point_value_old[i+EXPRESSION_VARIABLE_INDEX],&int_1);
			}

			// jacobian
			for(i = 0; i < term[t].n_variables; i ++)
			{
				v = term[t].variable[i];
				d = term[t].differential[i];

				expression_evaluate(n_points, point_term, term[t].jacobian[i], point_value, expression_work);
				for(p = 0; p < n_points; p ++) point_term[p] *= - element[e].W[p] * term[t].implicit;

				for(j = 0; j < n_basis[q]; j ++)
					for(p = 0; p < n_points; p ++)
						A[j][p] = element[e].P[x][j][p] * point_term[p];

				dgemm_(&trans[1],&trans[0],&n_basis[v],&n_basis[q],&n_points,
						&dbl_1,
						element[e].P[d][0],&n_points,
						A[0],&lda,
						&dbl_1,
						&local_value[0][local_element[q]][local_element[v]],&ldl);
			}

			// residual
			expression_evaluate(n_points, point_term, term[t].residual, point_value, expression_work);
			expression_evaluate(n_points, point_term_old, term[t].residual, point_value_old, expression_work);
			for(p = 0; p < n_points; p ++) point_term[p] = - element[e].W[p] *
				(term[t].implicit * point_term[p] + (1.0 - term[t].implicit) * point_term_old[p]);
			dgemv_(&trans[1],&n_points,&n_basis[q],
					&dbl_1,
					element[e].P[x][0],&n_points,
					point_term,&int_1,
					&dbl_1,
					&local_residual[0][local_element[q]],&int_1);
		}

		// add to the global system
		for(v = 0; v < n_variables; v ++)
		{
			for(i = 0; i < n_basis[v]; i ++)
			{
				sparse_add_to_row_values(system, element[e].unknown[v][i], local_value[0][local_element[v]+i]);
				residual[element[e].unknown[v][i]] += local_residual[0][local_element[v]+i];
			}
		}
	}

	// loop over the faces
	for(f = 0; f < n_faces; f ++)
	{
		// local system indices
		local_n = 0;
		for(v = 0; v < n_variables; v ++)
		{
			local_element[v] = local_n;
			local_n += n_basis[v];
		}
		for(b = 0; b < face[f].n_borders; b ++)
		{
			local_n = sum_n_basis;
			for(i = 0; i < face[f].border[b]->n_faces; i ++)
			{
				if(face[f].border[b]->face[i] == &face[f])
				{
					opposite[b] = i;
					for(v = 0; v < n_variables; v ++)
					{
						local_adjacent[b][v] = local_n;
						local_n += n_basis[v];
					}
				}
				else
				{
					for(j = 0; j < face[f].border[b]->face[i]->n_borders; j ++)
					{
						if(face[f].border[b]->face[i]->border[j] != face[f].border[b])
						{
							local_n += sum_n_basis;
						}
					}
				}
			}
		}

		// initialise the local system
		n = face[f].n_borders*sum_n_basis;
		dcopy_(&n,&dbl_0,&int_0,local_residual[0],&int_1);
		n = face[f].n_borders*sum_n_basis*ldl;
		dcopy_(&n,&dbl_0,&int_0,local_value[0][0],&int_1);

		// location and normal values at the intergration points
		for(i = 0; i < 2; i ++)
		{
			dcopy_(&n_gauss,face[f].X[i],&int_1,point_value[i+EXPRESSION_LOCATION_INDEX],&int_1);
			dcopy_(&n_gauss,face[f].X[i],&int_1,point_value_old[i+EXPRESSION_LOCATION_INDEX],&int_1);
			dcopy_(&n_gauss,&face[f].normal[i],&int_0,point_value[i+EXPRESSION_NORMAL_INDEX],&int_1);
			dcopy_(&n_gauss,&face[f].normal[i],&int_0,point_value_old[i+EXPRESSION_NORMAL_INDEX],&int_1);
		}

		// loop over the terms
		for(t = 0; t < n_terms; t ++)
		{
			if(term[t].type == 's') continue;

			// equation and direction
			q = term[t].equation;
			x = term[t].type == 'y';

			// variable values at the intergration points
			for(i = 0; i < term[t].n_variables; i ++)
			{
				v = term[t].variable[i];
				d = term[t].differential[i];

				for(j = 0; j < face[f].n_borders; j ++)
				{
					for(k = 0; k < n_basis[v]; k ++)
					{
						basis_value[k+j*n_basis[v]] = u[face[f].border[j]->unknown[v][k]];
						basis_value_old[k+j*n_basis[v]] = u_old[face[f].border[j]->unknown[v][k]];
					}
				}
				for(j = 0; j < face[f].n_boundaries[v]; j ++)
				{
					basis_value[j+face[f].n_borders*n_basis[v]] = face[f].boundary[v][j]->value;
					basis_value_old[j+face[f].n_borders*n_basis[v]] = face[f].boundary[v][j]->value;
				}

				if(term[t].method[i] == 'i' || d != powers_taylor[0][0] || face[f].n_borders < 2 || face[f].n_boundaries[v])
				{
					n = face[f].n_borders*n_basis[v] + face[f].n_boundaries[v];
					dgemv_(&trans[0],&n_gauss,&n,
							&dbl_1,
							face[f].Q[v][d][0],&n_gauss,
							basis_value,&int_1,
							&dbl_0,
							point_value[i+EXPRESSION_VARIABLE_INDEX],&int_1);
					dgemv_(&trans[0],&n_gauss,&n,
							&dbl_1,
							face[f].Q[v][d][0],&n_gauss,
							basis_value_old,&int_1,
							&dbl_0,
							point_value_old[i+EXPRESSION_VARIABLE_INDEX],&int_1);
				}
				else if(term[t].method[i] == 'a' || term[t].method[i] == 'd')
				{
					dcopy_(&n_gauss,&dbl_0,&int_0,point_value[i+EXPRESSION_VARIABLE_INDEX],&int_1);
					dcopy_(&n_gauss,&dbl_0,&int_0,point_value_old[i+EXPRESSION_VARIABLE_INDEX],&int_1);

					for(b = 0; b < face[f].n_borders; b ++)
					{
						alpha = term[t].method[i] == 'a' ? 0.5 : 0.5 * face[f].border[b]->orient[opposite[b]];
						dgemv_(&trans[0],&n_gauss,&n_basis[v],
								&alpha,
								face[f].border[b]->Q[opposite[b]][0],&n_gauss,
								&basis_value[b*n_basis[v]],&int_1,
								&dbl_1,
								point_value[i+EXPRESSION_VARIABLE_INDEX],&int_1);
						dgemv_(&trans[0],&n_gauss,&n_basis[v],
								&alpha,
								face[f].border[b]->Q[opposite[b]][0],&n_gauss,
								&basis_value_old[b*n_basis[v]],&int_1,
								&dbl_1,
								point_value_old[i+EXPRESSION_VARIABLE_INDEX],&int_1);
					}
				}
			}

			// jacobian
			for(i = 0; i < term[t].n_variables; i ++)
			{
				v = term[t].variable[i];
				d = term[t].differential[i];

				expression_evaluate(n_gauss, point_term, term[t].jacobian[i], point_value, expression_work);
				for(p = 0; p < n_gauss; p ++) point_term[p] *= face[f].normal[x] * term[t].implicit * face[f].W[p];

				for(b = 0; b < face[f].n_borders; b ++)
				{
					for(j = 0; j < n_basis[q]; j ++)
						for(p = 0; p < n_gauss; p ++)
							A[j][p] = face[f].border[b]->Q[opposite[b]][j][p] * point_term[p];

					if(term[t].method[i] == 'i' || d != powers_taylor[0][0] || face[f].n_borders < 2 || face[f].n_boundaries[v])
					{
						alpha = face[f].border[b]->orient[opposite[b]];
						for(j = 0; j < face[f].n_borders; j ++)
						{
							dgemm_(&trans[1],&trans[0],&n_basis[v],&n_basis[q],&n_gauss,
									&alpha,
									face[f].Q[v][d][j*n_basis[v]],&n_gauss,
									A[0],&lda,
									&dbl_1,
									&local_value[b][local_element[q]][j == b ? local_element[v] : local_adjacent[b][v]],&ldl);
						}
					}
					else if(term[t].method[i] == 'a' || term[t].method[i] == 'd')
					{
						alpha = 0.5 * face[f].border[b]->orient[opposite[b]];
						for(j = 0; j < face[f].n_borders; j ++)
						{
							dgemm_(&trans[1],&trans[0],&n_basis[v],&n_basis[q],&n_gauss,
									&alpha,
									face[f].border[j]->Q[opposite[j]][0],&n_gauss,
									A[0],&lda,
									&dbl_1,
									&local_value[b][local_element[q]][j == b ? local_element[v] : local_adjacent[b][v]],&ldl);
						}
					}
				}
			}

			// residual
			expression_evaluate(n_gauss, point_term, term[t].residual, point_value, expression_work);
			expression_evaluate(n_gauss, point_term_old, term[t].residual, point_value_old, expression_work);
			for(p = 0; p < n_gauss; p ++) point_term[p] = face[f].normal[x] * face[f].W[p] *
				(term[t].implicit * point_term[p] + (1.0 - term[t].implicit) * point_term_old[p]);
			for(b = 0; b < face[f].n_borders; b ++)
			{
				alpha = face[f].border[b]->orient[opposite[b]];
				dgemv_(&trans[1],&n_gauss,&n_basis[q],
						&alpha,
						face[f].border[b]->Q[opposite[b]][0],&n_gauss,
						point_term,&int_1,
						&dbl_1,
						&local_residual[b][local_element[q]],&int_1);
			}
		}

		// add to the global system
		for(b = 0; b < face[f].n_borders; b ++)
		{
			for(v = 0; v < n_variables; v ++)
			{
				for(i = 0; i < n_basis[v]; i ++)
				{
					sparse_add_to_row_values(system, face[f].border[b]->unknown[v][i], local_value[b][local_element[v]+i]);
					residual[face[f].border[b]->unknown[v][i]] += local_residual[b][local_element[v]+i];
				}
			}
		}
	}

	free(n_basis);

	destroy_matrix((void *)local_residual);
	destroy_tensor((void *)local_value);
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
