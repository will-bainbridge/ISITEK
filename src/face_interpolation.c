//////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "condition.h"
#include "face.h"
#include "linear.h"
#include "memory.h"
#include "numerics.h"
#include "solver.h"

//////////////////////////////////////////////////////////////////

#include "face.r"

//////////////////////////////////////////////////////////////////

static int n_variables, *order, max_order, *n_bases, max_n_bases, sum_n_bases;

static int n_interpolations, *interpolation_variable, *interpolation_differential;

static double **x, **y, ***x_adj, ***y_adj, **centre;

static double **R, **inv_R, **T;

static int *n_adj_bases, *n_int_bases;

static int *face_taylor;

static int ldp, lds, ldq, lda, ldb, ldf, ldd;
static int sizep;
static int incd, incq;
static double **P, **S, **Q, **A, **B, **F, ***D;

static int info, *pivot, int_1 = 1;
static char trans[2] = "NT";
static double d_one = 1.0, d_zero = 0.0;

//////////////////////////////////////////////////////////////////

int face_interpolation_start()
{
	// numbers
	n_variables = solver_n_variables();
	order = (int *)malloc(n_variables * sizeof(int));
	if(order == NULL) return FACE_MEMORY_ERROR;
	solver_variable_order(order);
	max_order = solver_variable_max_order();
	n_bases = (int *)malloc(n_variables * sizeof(int));
	if(n_bases == NULL) return FACE_MEMORY_ERROR;
	solver_variable_n_bases(n_bases);
	max_n_bases = solver_variable_max_n_bases();
	sum_n_bases = solver_variable_sum_n_bases();

	// quadrature
	int n_gauss = solver_n_gauss();
	int n_hammer = solver_n_hammer();

	// interplations
	n_interpolations = solver_n_interpolations();
	interpolation_variable = (int *)malloc(n_interpolations * sizeof(int));
	if(interpolation_variable == NULL) return FACE_MEMORY_ERROR;
	solver_interpolation_variable(interpolation_variable);
	interpolation_differential = (int *)malloc(n_interpolations * sizeof(int));
	if(interpolation_differential == NULL) return FACE_MEMORY_ERROR;
	solver_interpolation_differential(interpolation_differential);

	// locations
	y = matrix_double_new(NULL,2,n_gauss);
	x_adj = tensor_double_new(NULL,2,2,n_hammer*(ELEMENT_MAX_N_FACES-2));
	y_adj = tensor_double_new(NULL,2,2,n_hammer*(ELEMENT_MAX_N_FACES-2));
	centre = matrix_double_new(NULL,2,1);
	if(y == NULL || x_adj == NULL || y_adj == NULL || centre == NULL) return FACE_MEMORY_ERROR;

	// transformation
	R = matrix_double_new(NULL,2,2);
	inv_R = matrix_double_new(NULL,2,2);
	T = matrix_double_new(NULL,max_n_bases,max_n_bases);
	if(R == NULL || inv_R == NULL || T == NULL) return FACE_MEMORY_ERROR;

	// interpolation problem sizes
	int max_n_conditions = condition_max_n_variables();
	int max_n_adj_bases = 2*max_n_bases;
	int max_n_int_bases = 2*max_n_bases + max_n_conditions*n_gauss;
	n_adj_bases = (int *)malloc(n_variables * sizeof(int));
	if(n_adj_bases == NULL) return FACE_MEMORY_ERROR;
	n_int_bases = (int *)malloc(n_variables * sizeof(int));
	if(n_int_bases == NULL) return FACE_MEMORY_ERROR;

	// face basis indices
	face_taylor = (int *)malloc(max_n_int_bases * sizeof(int));
	if(face_taylor == NULL) return FACE_MEMORY_ERROR;

	// temporary matrices
	ldp = lds = ldq = 2*(ELEMENT_MAX_N_FACES-2)*n_hammer;
	lda = ldb = max_n_int_bases;
	ldf = ldd = n_gauss;
	sizep = max_n_adj_bases*ldp;
	incd = max_n_int_bases*n_gauss;
	P = matrix_double_new(NULL,max_n_adj_bases,ldp);
	S = matrix_double_new(NULL,max_n_adj_bases,lds);
	Q = matrix_double_new(NULL,max_n_int_bases,ldp);
	A = matrix_double_new(NULL,max_n_int_bases,max_n_int_bases);
	B = matrix_double_new(NULL,max_n_int_bases,max_n_int_bases);
	F = matrix_double_new(NULL,max_n_int_bases,n_gauss);
	D = tensor_double_new(NULL,max_n_bases,max_n_int_bases,n_gauss);
	if(P == NULL || S == NULL || Q == NULL || A == NULL || B == NULL || F == NULL || D == NULL) return FACE_MEMORY_ERROR;

	return ELEMENT_SUCCESS;
}

//////////////////////////////////////////////////////////////////

void face_interpolation_end()
{
	free(order);
	free(n_bases);
	free(interpolation_variable);
	free(interpolation_differential);
	matrix_free((void *)y);
	tensor_free((void *)x_adj);
	tensor_free((void *)y_adj);
	matrix_free((void *)centre);
	matrix_free((void *)R);
	matrix_free((void *)inv_R);
	matrix_free((void *)T);
	free(n_adj_bases);
	free(n_int_bases);
	free(face_taylor);
	matrix_free((void *)P);
	matrix_free((void *)S);
	matrix_free((void *)Q);
	matrix_free((void *)A);
	matrix_free((void *)B);
	matrix_free((void *)F);
	tensor_free((void *)D);
}

//////////////////////////////////////////////////////////////////

int face_interpolation_calculate(FACE face)
{
	int a, c, i, j, q, v;

	CONDITION condition = face->boundary ? boundary_condition(face->boundary) : NULL;
	int n_conditions = condition ? condition_n_variables(condition) : 0;

	for(v = 0; v < n_variables; v ++) n_adj_bases[v] = n_int_bases[v] = face->n_borders*n_bases[v];
	for(c = 0; c < n_conditions; c ++) n_int_bases[condition_variable(condition,c)] += face->n_quadrature;

	int sum_n_adj_bases = 0, sum_n_int_bases = 0;
	for(i = 0; i < n_interpolations; i ++)
	{
		sum_n_adj_bases += n_adj_bases[interpolation_variable[i]];
		sum_n_int_bases += n_int_bases[interpolation_variable[i]];
	}

	// allocate
	face->Q = (double ***)malloc(n_interpolations * sizeof(double **));
	if(face->Q == NULL) return FACE_MEMORY_ERROR;
	face->Q[0] = (double **)malloc(sum_n_int_bases * sizeof(double *));
	if(face->Q[0] == NULL) return FACE_MEMORY_ERROR;
	face->Q[0][0] = (double *)malloc(sum_n_int_bases * face->n_quadrature * sizeof(double));
	if(face->Q[0][0] == NULL) return FACE_MEMORY_ERROR;
	for(i = 1; i < n_interpolations; i ++) face->Q[i] = face->Q[i-1] + n_int_bases[interpolation_variable[i-1]];
	for(i = 1; i < n_interpolations; i ++) face->Q[i][0] = face->Q[i-1][0] + n_int_bases[interpolation_variable[i-1]] * face->n_quadrature;
	for(i = 0; i < n_interpolations; i ++) for(j = 1; j < n_int_bases[interpolation_variable[i]]; j ++) face->Q[i][j] = face->Q[i][j-1] + face->n_quadrature;

	// rotation to face coordinates
	R[0][0] = + face->normal[0]; R[0][1] = + face->normal[1];
	R[1][0] = - face->normal[1]; R[1][1] = + face->normal[0];
	double det_R = R[0][0]*R[1][1] - R[0][1]*R[1][0];
	inv_R[0][0] = + R[1][1]/det_R; inv_R[0][1] = - R[0][1]/det_R;
	inv_R[1][0] = - R[1][0]/det_R; inv_R[1][1] = + R[0][0]/det_R;
	numerics_transformation_matrix(max_order,T,inv_R);

	// face integration locations
	x = face->X;
	for(q = 0; q < face->n_quadrature; q ++)
	{
		for(i = 0; i < 2; i ++)
		{
			y[i][q] = face->centre[i];
			for(j = 0; j < 2; j ++) y[i][q] += R[i][j]*(x[j][q] - face->centre[j]);
		}
	}

	// numbers of element integration locations
	int n_points[2], sum_n_points[3];
	for(a = 0; a < face->n_borders; a ++) n_points[a] = element_n_quadrature(face->border[a]);
	sum_n_points[0] = 0;
	for(a = 0; a < face->n_borders; a ++) sum_n_points[a+1] = sum_n_points[a] + n_points[a];

	// adjacent element integration locations
	for(a = 0; a < face->n_borders; a ++)
	{
		element_quadrature_x(face->border[a],x_adj[a]);
		for(q = 0; q < n_points[a]; q ++)
		{
			for(i = 0; i < 2; i ++)
			{
				y_adj[a][i][q] = face->centre[i];
				for(j = 0; j < 2; j ++) y_adj[a][i][q] += R[i][j]*(x_adj[a][j][q] - face->centre[j]);
			}
		}
	}

	/*// for all variables
	for(v = 0; v < n_variables; v ++)
	{
		if(!max_update && !update[v]) continue;

		n_adj = face[f].n_borders;
		adj = face[f].border;
		n_bnd = face[f].n_boundaries[v];
		bnd = face[f].boundary[v];

		n_adj_bases = n_adj*n_basis[v];
		n_int_terms = n_adj_bases + n_bnd*n_gauss;

		// face basis indices
		n_int_bases = 0;
		for(i = 0; i < n_adj*variable_order[v] + n_bnd; i ++)
			for(j = 0; j < n_adj*variable_order[v] + n_bnd; j ++)
				if(i + n_adj*j < n_adj*variable_order[v] + n_bnd && j < n_gauss)
					face_taylor[n_int_bases ++] = powers_taylor[i][j];
		exit_if_false(n_int_bases == n_adj_bases + n_bnd*n_gauss,"mismatched number of interpolation unknowns");

		// element bases at the integration locations
		for(i = 0; i < n_adj*n_basis[v]; i ++) for(j = 0; j < sum_n_points[n_adj]; j ++) P[i][j] = 0.0;
		for(a = 0; a < n_adj; a ++)
			for(i = 0; i < n_basis[v]; i ++)
				basis(n_points[a],&P[i+a*n_basis[v]][sum_n_points[a]],x_adj[a],adj[a]->centre,adj[a]->size,i,none);

		// face bases at the integration locations
		for(a = 0; a < n_adj; a ++)
			for(i = 0; i < n_int_bases; i ++)
				basis(n_points[a],&Q[i][sum_n_points[a]],y_adj[a],face[f].centre,face[f].size,face_taylor[i],none);

		// centre of face in form which can be passed to basis
		for(i = 0; i < 2; i ++) centre[i][0] = face[f].centre[i];

		// integration matrix
		dcopy_(&sizep,P[0],&i_one,S[0],&i_one);
		for(a = 0; a < n_adj; a ++)
			for(i = 0; i < n_points[a]; i ++)
				dscal_(&n_basis[v],&adj[a]->W[i],&S[a*n_basis[v]][i+sum_n_points[a]],&lds);

		// weak interpolation system
		dgemm_(&trans[1],&trans[0],&n_adj_bases,&n_int_bases,&sum_n_points[n_adj],&d_one,S[0],&lds,Q[0],&ldq,&d_zero,A[0],&lda);

		// weak interpolation rhs
		dgemm_(&trans[1],&trans[0],&n_adj_bases,&n_adj_bases,&sum_n_points[n_adj],&d_one,S[0],&lds,P[0],&ldp,&d_zero,B[0],&ldb);

		// boundary conditions
		for(b = 0; b < n_bnd; b ++)
		{
			for(i = 0; i < n_int_bases; i ++)
				basis(n_gauss,&A[i][n_adj_bases],y,face[f].centre,face[f].size,face_taylor[i],taylor_powers[bnd[b]->condition]);

			for(i = 0; i < n_gauss; i ++)
				for(j = 0; j < n_int_terms; j ++)
					B[j][i+b*n_gauss+n_adj_bases] = B[i+b*n_gauss+n_adj_bases][j] = (i+b*n_gauss+n_adj_bases) == j;
		}

		//if(n_bnd)
		//{
		//	printf("\nface %i variable %i\n\n",f,v);
		//
		//	for(i = 0; i < n_gauss; i ++) { printf("(%lf %lf) ",face[f].X[0][i],face[f].X[1][i]); } printf("\n\n");
		//
		//	for(i = 0; i < n_int_bases; i ++) { printf("(%i,%i) ",taylor_powers[face_taylor[i]][0],taylor_powers[face_taylor[i]][1]); } printf("\n\n");
		//
		//	for(i = 0; i < n_int_bases; i ++) { for(j = 0; j < n_int_bases; j ++) { printf("%+.1e ",A[j][i]); } printf("\n"); }
		//	for(i = 0; i < n_int_bases*9-1; i ++) { printf("-"); } printf("\n");
		//	for(i = 0; i < n_int_bases; i ++) { for(j = 0; j < n_int_terms; j ++) { printf("%+.1e ",B[j][i]); } printf("\n"); } printf("\n");
		//	
		//	getchar();
		//}

		// solve interpolation problem
		dgesv_(&n_int_bases,&n_int_terms,A[0],&lda,pivot,B[0],&ldb,&info);

		// interpolate values to the face integration locations
		for(i = 0; i < n_basis[v]; i ++)
		{
			for(j = 0; j < n_int_bases; j ++) basis(n_gauss,F[j],y,face[f].centre,face[f].size,face_taylor[j],taylor_powers[i]);
			dgemm_(&trans[0],&trans[0],&n_gauss,&n_int_terms,&n_int_bases,&d_one,F[0],&ldf,B[0],&ldb,&d_zero,D[i][0],&ldd);
		}

		// transform from face to cartesian coordinates
		incq = n_int_terms*n_gauss;
		for(i = 0; i < n_int_terms; i ++)
			for(j = 0; j < n_gauss; j ++)
				dgemv_(&trans[1],&n_basis[v],&n_basis[v],&d_one,T[0],&max_n_basis,&D[0][i][j],&incd,&d_zero,&face[f].Q[v][0][i][j],&incq);

		//if(n_bnd)
		//{
		//	for(j = 0; j < n_gauss; j ++) {
		//		for(i = 0; i < n_int_terms; i ++) {
		//			printf("%+.1e ",D[0][i][j]);
		//		} printf("\n");
		//	} printf("\n");
		//	getchar();
		//}
	}*/

	return FACE_SUCCESS;
}

//////////////////////////////////////////////////////////////////

/*int update_face_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face, int n_boundaries_old, struct BOUNDARY *boundary_old)
{
	int a, b, f, g, h, i, j, v;

	// lists of old face boundaries
	int **face_n_boundaries_old = allocate_integer_matrix(NULL,n_faces,n_variables);
	exit_if_false(face_n_boundaries_old != NULL,"allocating face_n_boundaries_old");
	for(i = 0; i < n_faces; i ++) for(j = 0; j < n_variables; j ++) face_n_boundaries_old[i][j] = 0;

	int ***face_boundary_old = allocate_integer_tensor(NULL,n_faces,n_variables,MAX_FACE_N_BOUNDARIES);
	exit_if_false(face_boundary_old != NULL,"allocating face_boundary_old");

	for(b = 0; b < n_boundaries_old; b ++)
	{
		for(i = 0; i < boundary_old[b].n_faces; i ++)
		{
			v = boundary_old[b].variable;
			f = boundary_old[b].face[i] - &face[0];
			face_boundary_old[f][v][face_n_boundaries_old[f][v]++] = b;
		}
	}

	// order and numbers of integration points and bases
	int max_variable_order = 0, max_variable_order_old = 0;
	for(v = 0; v < n_variables; v ++) max_variable_order = MAX(max_variable_order,variable_order[v]);
	for(v = 0; v < n_variables_old; v ++) max_variable_order_old = MAX(max_variable_order_old,variable_order_old[v]);

	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order);
	int n_points[MAX_FACE_N_BORDERS], sum_n_points[MAX_FACE_N_BORDERS+1];

	int *n_basis = (int *)malloc(n_variables * sizeof(int)), max_n_basis = ORDER_TO_N_BASIS(max_variable_order);
	exit_if_false(n_basis != NULL,"allocating n_basis");
	for(v = 0; v < n_variables; v ++) n_basis[v] = ORDER_TO_N_BASIS(variable_order[v]);

	// truth values for updating the interpolation on a face
	int *update = (int *)malloc(n_variables * sizeof(int)), max_update = max_variable_order != max_variable_order_old, any_update;
	exit_if_false(update != NULL,"allocating update");

	// transformation
	double **R = allocate_double_matrix(NULL,2,2), det_R, **inv_R = allocate_double_matrix(NULL,2,2);
	exit_if_false(R != NULL,"allocating R");
	exit_if_false(inv_R != NULL,"allocating inv_R");
	double **T = allocate_double_matrix(NULL,max_n_basis,max_n_basis);
	exit_if_false(T != NULL,"allocating T");

	// centre
	double **centre = allocate_double_matrix(NULL,2,1);
	exit_if_false(centre != NULL,"allocating centre");

	// zero-differential
	int none[2] = {0,0};

	// integration locations in cartesian (x) and tranformed (y) coordinates
	double **x, **y;
	double **x_adj[MAX_FACE_N_BORDERS], ***y_adj;
	exit_if_false(y = allocate_double_matrix(NULL,2,n_gauss),"allocating y");
	exit_if_false(y_adj = allocate_double_tensor(NULL,MAX_FACE_N_BORDERS,2,n_hammer*(MAX_ELEMENT_N_FACES-2)),"allocating y_adj");

	// adjacent elements and boundaries to the face
	int n_adj, n_bnd;
	struct ELEMENT **adj;
	struct BOUNDARY **bnd;

	// interpolation problem sizes
	int n_adj_bases = MAX_FACE_N_BORDERS*max_n_basis; // number of adjacent bases
	int n_int_terms = MAX_FACE_N_BORDERS*max_n_basis + MAX_FACE_N_BOUNDARIES; // number of terms from which to interpolate
	int n_int_bases = MAX_FACE_N_BORDERS*max_n_basis + MAX_FACE_N_BOUNDARIES*max_variable_order; // number of interpolant bases

	// face basis taylor indices
	int *face_taylor = (int *)malloc(n_int_bases * sizeof(int));
	exit_if_false(face_taylor != NULL,"allocating face_taylor");

	// temporary matrices
	int ldp, lds, ldq;
	ldp = lds = ldq = MAX_FACE_N_BORDERS*(MAX_ELEMENT_N_FACES-2)*n_hammer;
	int sizep = n_adj_bases*ldp;
	double **P = allocate_double_matrix(NULL,n_adj_bases,ldp);
	double **S = allocate_double_matrix(NULL,n_adj_bases,lds);
	double **Q = allocate_double_matrix(NULL,n_int_bases,ldp);

	int lda, ldb;
	lda = ldb = n_int_bases;
	double **A = allocate_double_matrix(NULL,n_int_bases,n_int_bases);
	double **B = allocate_double_matrix(NULL,n_int_bases,n_int_bases);

	int ldf, ldd;
	ldf = ldd = n_gauss;
	double **F = allocate_double_matrix(NULL,n_int_bases,n_gauss);
	double ***D = allocate_double_tensor(NULL,max_n_basis,n_int_terms,n_gauss);

	int incd = n_int_terms*n_gauss, incq;

	exit_if_false(P != NULL && S != NULL && Q != NULL && A != NULL && B != NULL && F != NULL && D != NULL,"allocating working matrices");

	// blas 
	char trans[2] = "NT";
	int i_one = 1;
	double d_one = 1.0, d_zero = 0.0;

	// lapack
	int info, *pivot = (int *)malloc(n_int_bases * sizeof(int));
	exit_if_false(pivot != NULL,"allocating pivot");

	// number of updated interpolations
	int updated = 0;

	for(f = 0; f < n_faces; f ++)
	{
		// see if face interpolation needs updating
		any_update = 0;
		for(v = 0; v < n_variables; v ++)
		{
			update[v] = 0;
			if(v >= n_variables_old) update[v] = 1;
			else if(variable_order[v] != variable_order_old[v]) update[v] = 1;
			else if(face[f].n_boundaries[v] != face_n_boundaries_old[f][v]) update[v] = 1;
			else for(i = 0; i < face[f].n_boundaries[v]; i ++)
				if(face[f].boundary[v][i]->condition != boundary_old[face_boundary_old[f][v][i]].condition)
					update[v] = 1;
			any_update += update[v];
		}
		if(!any_update) continue;

		// allocate matrices
		exit_if_false(face[f].Q = allocate_face_q(&face[f],n_variables,n_basis,n_gauss),"allocating face Q");

		// rotation to face coordinates
		R[0][0] = + face[f].normal[0]; R[0][1] = + face[f].normal[1];
		R[1][0] = - face[f].normal[1]; R[1][1] = + face[f].normal[0];
		det_R = R[0][0]*R[1][1] - R[0][1]*R[1][0];
		inv_R[0][0] = + R[1][1]/det_R; inv_R[0][1] = - R[0][1]/det_R;
		inv_R[1][0] = - R[1][0]/det_R; inv_R[1][1] = + R[0][0]/det_R;
		transformation_matrix(max_variable_order,T,inv_R);

		// face integration locations
		x = face[f].X;
		for(g = 0; g < n_gauss; g ++)
		{
			for(i = 0; i < 2; i ++)
			{
				y[i][g] = face[f].centre[i];
				for(j = 0; j < 2; j ++) y[i][g] += R[i][j]*(x[j][g] - face[f].centre[j]);
			}
		}

		// numbers of element integration locations
		for(a = 0; a < face[f].n_borders; a ++) n_points[a] = n_hammer*(face[f].border[a]->n_faces-2);
		sum_n_points[0] = 0;
		for(a = 0; a < face[f].n_borders; a ++) sum_n_points[a+1] = sum_n_points[a] + n_points[a];

		// adjacent element integration locations
		for(a = 0; a < face[f].n_borders; a ++)
		{
			x_adj[a] = face[f].border[a]->X;
			for(h = 0; h < n_points[a]; h ++)
			{
				for(i = 0; i < 2; i ++)
				{
					y_adj[a][i][h] = face[f].centre[i];
					for(j = 0; j < 2; j ++) y_adj[a][i][h] += R[i][j]*(x_adj[a][j][h] - face[f].centre[j]);
				}
			}
		}

		// for all variables
		for(v = 0; v < n_variables; v ++)
		{
			if(!max_update && !update[v]) continue;

			n_adj = face[f].n_borders;
			adj = face[f].border;
			n_bnd = face[f].n_boundaries[v];
			bnd = face[f].boundary[v];

			n_adj_bases = n_adj*n_basis[v];
			n_int_terms = n_adj_bases + n_bnd*n_gauss;

			// face basis indices
			n_int_bases = 0;
			for(i = 0; i < n_adj*variable_order[v] + n_bnd; i ++)
				for(j = 0; j < n_adj*variable_order[v] + n_bnd; j ++)
					if(i + n_adj*j < n_adj*variable_order[v] + n_bnd && j < n_gauss)
						face_taylor[n_int_bases ++] = powers_taylor[i][j];
			exit_if_false(n_int_bases == n_adj_bases + n_bnd*n_gauss,"mismatched number of interpolation unknowns");

			// element bases at the integration locations
			for(i = 0; i < n_adj*n_basis[v]; i ++) for(j = 0; j < sum_n_points[n_adj]; j ++) P[i][j] = 0.0;
			for(a = 0; a < n_adj; a ++)
				for(i = 0; i < n_basis[v]; i ++)
					basis(n_points[a],&P[i+a*n_basis[v]][sum_n_points[a]],x_adj[a],adj[a]->centre,adj[a]->size,i,none);

			// face bases at the integration locations
			for(a = 0; a < n_adj; a ++)
				for(i = 0; i < n_int_bases; i ++)
					basis(n_points[a],&Q[i][sum_n_points[a]],y_adj[a],face[f].centre,face[f].size,face_taylor[i],none);

			// centre of face in form which can be passed to basis
			for(i = 0; i < 2; i ++) centre[i][0] = face[f].centre[i];

			// integration matrix
			dcopy_(&sizep,P[0],&i_one,S[0],&i_one);
			for(a = 0; a < n_adj; a ++)
				for(i = 0; i < n_points[a]; i ++)
					dscal_(&n_basis[v],&adj[a]->W[i],&S[a*n_basis[v]][i+sum_n_points[a]],&lds);

			// weak interpolation system
			dgemm_(&trans[1],&trans[0],&n_adj_bases,&n_int_bases,&sum_n_points[n_adj],&d_one,S[0],&lds,Q[0],&ldq,&d_zero,A[0],&lda);

			// weak interpolation rhs
			dgemm_(&trans[1],&trans[0],&n_adj_bases,&n_adj_bases,&sum_n_points[n_adj],&d_one,S[0],&lds,P[0],&ldp,&d_zero,B[0],&ldb);

			// boundary conditions
			for(b = 0; b < n_bnd; b ++)
			{
				for(i = 0; i < n_int_bases; i ++)
					basis(n_gauss,&A[i][n_adj_bases],y,face[f].centre,face[f].size,face_taylor[i],taylor_powers[bnd[b]->condition]);

				for(i = 0; i < n_gauss; i ++)
					for(j = 0; j < n_int_terms; j ++)
						B[j][i+b*n_gauss+n_adj_bases] = B[i+b*n_gauss+n_adj_bases][j] = (i+b*n_gauss+n_adj_bases) == j;
			}

			//if(n_bnd)
			//{
			//	printf("\nface %i variable %i\n\n",f,v);
			//
			//	for(i = 0; i < n_gauss; i ++) { printf("(%lf %lf) ",face[f].X[0][i],face[f].X[1][i]); } printf("\n\n");
			//
			//	for(i = 0; i < n_int_bases; i ++) { printf("(%i,%i) ",taylor_powers[face_taylor[i]][0],taylor_powers[face_taylor[i]][1]); } printf("\n\n");
			//
			//	for(i = 0; i < n_int_bases; i ++) { for(j = 0; j < n_int_bases; j ++) { printf("%+.1e ",A[j][i]); } printf("\n"); }
			//	for(i = 0; i < n_int_bases*9-1; i ++) { printf("-"); } printf("\n");
			//	for(i = 0; i < n_int_bases; i ++) { for(j = 0; j < n_int_terms; j ++) { printf("%+.1e ",B[j][i]); } printf("\n"); } printf("\n");
			//	
			//	getchar();
			//}

			// solve interpolation problem
			dgesv_(&n_int_bases,&n_int_terms,A[0],&lda,pivot,B[0],&ldb,&info);

			// interpolate values to the face integration locations
			for(i = 0; i < n_basis[v]; i ++)
			{
				for(j = 0; j < n_int_bases; j ++) basis(n_gauss,F[j],y,face[f].centre,face[f].size,face_taylor[j],taylor_powers[i]);
				dgemm_(&trans[0],&trans[0],&n_gauss,&n_int_terms,&n_int_bases,&d_one,F[0],&ldf,B[0],&ldb,&d_zero,D[i][0],&ldd);
			}

			// transform from face to cartesian coordinates
			incq = n_int_terms*n_gauss;
			for(i = 0; i < n_int_terms; i ++)
				for(j = 0; j < n_gauss; j ++)
					dgemv_(&trans[1],&n_basis[v],&n_basis[v],&d_one,T[0],&max_n_basis,&D[0][i][j],&incd,&d_zero,&face[f].Q[v][0][i][j],&incq);

			//if(n_bnd)
			//{
			//	for(j = 0; j < n_gauss; j ++) {
			//		for(i = 0; i < n_int_terms; i ++) {
			//			printf("%+.1e ",D[0][i][j]);
			//		} printf("\n");
			//	} printf("\n");
			//	getchar();
			//}

			updated ++;
		}
	}

	// clean up
	destroy_matrix((void *)face_n_boundaries_old);
	destroy_tensor((void *)face_boundary_old);

	free(n_basis);
	free(update);

	destroy_matrix((void *)R);
	destroy_matrix((void *)inv_R);
	destroy_matrix((void *)T);

	destroy_matrix((void *)centre);
	destroy_matrix((void *)y);
	destroy_tensor((void *)y_adj);

	free(face_taylor);

	destroy_matrix((void *)P);
	destroy_matrix((void *)Q);
	destroy_matrix((void *)S);
	destroy_matrix((void *)A);
	destroy_matrix((void *)B);
	destroy_matrix((void *)F);
	destroy_tensor((void *)D);

	free(pivot);

	return updated;
}*/

//////////////////////////////////////////////////////////////////
