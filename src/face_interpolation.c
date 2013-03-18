//////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "condition.h"
#include "face.h"
#include "info.h"
#include "linear.h"
#include "memory.h"
#include "numerics.h"
#include "solver.h"

//////////////////////////////////////////////////////////////////

#include "face.r"

//////////////////////////////////////////////////////////////////

static int n_variables, max_n_bases, sum_n_bases;

static int n_gauss, n_hammer;

static int n_interpolations, *interpolation_variable, *interpolation_differential;

static int *n_constraints, **constraint_temporary, **constraint_differential;

static double **x, **y, ***y_adj;

static double **R, **inv_R, **T;

static int *n_adj_bases, *n_int_bases;

static int *face_taylor;

static int ldp, lds, ldq, lda, ldb, ldf, ldd;
static int sizep, sizef;
static double **P, **S, **Q, **A, **B, **F, ***D;

static int info, *pivot, int_1 = 1;
static char trans[2] = "NT";
static double double_1 = 1.0, double_0 = 0.0;

static int taylor_power[2];

//////////////////////////////////////////////////////////////////

int face_interpolation_start()
{
	// numbers
	n_variables = solver_n_variables();
	max_n_bases = solver_variable_max_n_bases();
	sum_n_bases = solver_variable_sum_n_bases();

	// quadrature
	n_gauss = solver_n_gauss();
	n_hammer = solver_n_hammer();

	// interplations
	n_interpolations = solver_n_interpolations();
	interpolation_variable = (int *)malloc(n_interpolations * sizeof(int));
	if(interpolation_variable == NULL) return FACE_MEMORY_ERROR;
	solver_interpolation_variable(interpolation_variable);
	interpolation_differential = (int *)malloc(n_interpolations * sizeof(int));
	if(interpolation_differential == NULL) return FACE_MEMORY_ERROR;
	solver_interpolation_differential(interpolation_differential);

	// constraints
	n_constraints = (int *)malloc(n_variables * sizeof(int));
	constraint_temporary = matrix_integer_new(NULL,2,condition_max_n_variables());
	constraint_differential = matrix_integer_new(NULL,n_variables,condition_max_n_variables());
	if(n_constraints == NULL || constraint_temporary == NULL || constraint_differential == NULL) return FACE_MEMORY_ERROR;

	// locations
	y = matrix_double_new(NULL,2,n_gauss);
	y_adj = tensor_double_new(NULL,2,2,n_hammer*(ELEMENT_MAX_N_FACES-2));
	if(y == NULL || y_adj == NULL) return FACE_MEMORY_ERROR;

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
	sizef = max_n_int_bases*ldf;
	P = matrix_double_new(NULL,max_n_adj_bases,ldp);
	S = matrix_double_new(NULL,max_n_adj_bases,lds);
	Q = matrix_double_new(NULL,max_n_int_bases,ldp);
	A = matrix_double_new(NULL,max_n_int_bases,max_n_int_bases);
	B = matrix_double_new(NULL,max_n_int_bases,max_n_int_bases);
	F = matrix_double_new(NULL,max_n_int_bases,n_gauss);
	D = tensor_double_new(NULL,max_n_bases,max_n_int_bases,n_gauss);
	if(P == NULL || S == NULL || Q == NULL || A == NULL || B == NULL || F == NULL || D == NULL) return FACE_MEMORY_ERROR;

	// lapack
	pivot = (int *)malloc(max_n_int_bases * sizeof(int));
	if(pivot == NULL) return FACE_MEMORY_ERROR;

	return ELEMENT_SUCCESS;
}

//////////////////////////////////////////////////////////////////

void face_interpolation_end()
{
	free(interpolation_variable);
	free(interpolation_differential);
	free(n_constraints);
	matrix_free((void *)constraint_temporary);
	matrix_free((void *)constraint_differential);
	matrix_free((void *)y);
	tensor_free((void *)y_adj);
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
	free(pivot);
}

//////////////////////////////////////////////////////////////////

int face_interpolation_calculate(FACE face)
{
	int a, c, i, j, n, q, v;

	// boundary condition constraints
	CONDITION condition = face->boundary ? boundary_condition(face->boundary) : NULL;
	for(v = 0; v < n_variables; v ++) n_constraints[v] = 0;
	if(condition)
	{
		condition_variable(condition,constraint_temporary[0]);
		condition_differential(condition,constraint_temporary[1]);
		for(c = 0; c < condition_n_variables(condition); c ++)
			constraint_differential[constraint_temporary[0][c]][n_constraints[constraint_temporary[0][c]]++] = constraint_temporary[1][c];
	}

	// numbers of bases
	int sum_n_adj_bases = 0, sum_n_int_bases = 0;
	for(v = 0; v < n_variables; v ++)
	{
		n_adj_bases[v] = face->n_borders*solver_variable_n_bases()[v];
		n_int_bases[v] = n_adj_bases[v] + n_constraints[v]*face->n_quadrature;
	}
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
	numerics_transformation_matrix(solver_variable_max_order(),T,inv_R);

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

	// adjacent element geometry
	for(a = 0; a < face->n_borders; a ++)
	{
		for(q = 0; q < n_points[a]; q ++)
		{
			for(i = 0; i < 2; i ++)
			{
				y_adj[a][i][q] = face->centre[i];
				for(j = 0; j < 2; j ++) y_adj[a][i][q] += R[i][j]*(element_quadrature_x(face->border[a])[j][q] - face->centre[j]);
			}
		}
	}

	// for all variables
	for(v = 0; v < n_variables; v ++)
	{
		// face basis indices
		n = 0;
		for(i = 0; i < face->n_borders*solver_variable_order()[v] + n_constraints[v]; i ++)
			for(j = 0; j < face->n_borders*solver_variable_order()[v] + n_constraints[v]; j ++)
				if(i + face->n_borders*j < face->n_borders*solver_variable_order()[v] + n_constraints[v] && j < n_gauss)
					face_taylor[n ++] = numerics_power_taylor(i,j);
		exit_if_false(n == n_int_bases[v],"mismatched number of interpolation bases");

		// no differential
		for(i = 0; i < 2; i ++) taylor_power[i] = 0;

		// element bases at the integration locations
		for(i = 0; i < face->n_borders*solver_variable_n_bases()[v]; i ++) for(j = 0; j < sum_n_points[face->n_borders]; j ++) P[i][j] = 0.0;
		for(a = 0; a < face->n_borders; a ++)
			for(i = 0; i < solver_variable_n_bases()[v]; i ++)
				numerics_basis(n_points[a],&P[i+a*solver_variable_n_bases()[v]][sum_n_points[a]],element_quadrature_x(face->border[a]),element_centre(face->border[a]),element_size(face->border[a]),i,taylor_power);

		// face bases at the integration locations
		for(a = 0; a < face->n_borders; a ++)
			for(i = 0; i < n_int_bases[v]; i ++)
				numerics_basis(n_points[a],&Q[i][sum_n_points[a]],y_adj[a],face->centre,0.5*face->size,face_taylor[i],taylor_power);

		// integration matrix
		dcopy_(&sizep,P[0],&int_1,S[0],&int_1);
		for(a = 0; a < face->n_borders; a ++)
			for(i = 0; i < n_points[a]; i ++)
				dscal_(&solver_variable_n_bases()[v],&element_quadrature_w(face->border[a])[i],&S[a*solver_variable_n_bases()[v]][i+sum_n_points[a]],&lds);

		// weak interpolation system
		dgemm_(&trans[1],&trans[0],&n_adj_bases[v],&n_int_bases[v],&sum_n_points[face->n_borders],&double_1,S[0],&lds,Q[0],&ldq,&double_0,A[0],&lda);

		// weak interpolation rhs
		dgemm_(&trans[1],&trans[0],&n_adj_bases[v],&n_adj_bases[v],&sum_n_points[face->n_borders],&double_1,S[0],&lds,P[0],&ldp,&double_0,B[0],&ldb);

		// boundary conditions
		for(c = 0; c < n_constraints[v]; c ++)
		{
			for(i = 0; i < 2; i ++) taylor_power[i] = numerics_taylor_power(constraint_differential[v][c],i);
			for(i = 0; i < n_int_bases[v]; i ++)
				numerics_basis(n_gauss,&A[i][n_adj_bases[v]],y,face->centre,0.5*face->size,face_taylor[i],taylor_power);

			for(i = 0; i < n_gauss; i ++)
				for(j = 0; j < n_int_bases[v]; j ++)
					B[j][i+c*n_gauss+n_adj_bases[v]] = B[i+c*n_gauss+n_adj_bases[v]][j] = (i+c*n_gauss+n_adj_bases[v]) == j;
		}

		//if(condition)
		//{
		//	printf("\nface %i variable %i\n",face->index,v);
		//	for(i = 0; i < n_gauss; i ++) { printf("(%lf %lf)",face->X[0][i],face->X[1][i]); } printf("\n");
		//	for(i = 0; i < n_int_bases[v]; i ++) { printf("(%i,%i)",numerics_taylor_power(face_taylor[i],0),numerics_taylor_power(face_taylor[i],1)); } printf("\n");
		//	for(i = 0; i < n_int_bases[v]*9-1; i ++) { printf("-"); } printf("\n");
		//	for(i = 0; i < n_int_bases[v]; i ++) { for(j = 0; j < n_int_bases[v]; j ++) { printf("%+.1e ",A[j][i]); } printf("\n"); }
		//	for(i = 0; i < n_int_bases[v]*9-1; i ++) { printf("-"); } printf("\n");
		//	for(i = 0; i < n_int_bases[v]; i ++) { for(j = 0; j < n_int_bases[v]; j ++) { printf("%+.1e ",B[j][i]); } printf("\n"); } printf("\n");
		//	getchar();
		//}

		// solve interpolation problem
		dgesv_(&n_int_bases[v],&n_int_bases[v],A[0],&lda,pivot,B[0],&ldb,&info);

		// interpolate values to the face integration locations
		for(i = 0; i < solver_variable_n_bases()[v]; i ++)
		{
			for(j = 0; j < 2; j ++) taylor_power[j] = numerics_taylor_power(i,j);
			for(j = 0; j < n_int_bases[v]; j ++) numerics_basis(n_gauss,F[j],y,face->centre,0.5*face->size,face_taylor[j],taylor_power);
			dgemm_(&trans[0],&trans[0],&n_gauss,&n_int_bases[v],&n_int_bases[v],&double_1,F[0],&ldf,B[0],&ldb,&double_0,D[i][0],&ldd);
		}

		/*// transform from face to cartesian coordinates
		for(i = 0; i < n_interpolations; i ++)
		{
			if(interpolation_variable[i] == v)
			{
				for(j = 0; j < n_int_bases[v]; j ++)
				{
					int k ,l;
					for(k = 0; k < face->n_quadrature; k ++)
					{
						face->Q[i][j][k] = 0;
						for(l = 0; l < solver_variable_n_bases()[v]; l ++) face->Q[i][j][k] += T[interpolation_differential[i]][l]*D[l][j][k];
					}
				}
			}
		}*/

		// transform from face to cartesian coordinates
		n = n_int_bases[v]*face->n_quadrature;
		for(i = 0; i < n_interpolations; i ++)
			if(interpolation_variable[i] == v)
				dgemv_(&trans[0],&n,&solver_variable_n_bases()[v],&double_1,D[0][0],&sizef,T[interpolation_differential[i]],&int_1,&double_0,face->Q[i][0],&int_1);

		//if(condition)
		//{
		//	for(i = 0; i < solver_variable_n_bases()[v]; i ++) {
		//		for(j = 0; j < solver_variable_n_bases()[v]; j ++) {
		//			printf("%+.1e ",T[i][j]);
		//		} printf("\n");
		//	} printf("\n");
		//	for(i = 0; i < solver_variable_n_bases()[v]; i ++) {
		//		for(j = 0; j < n_gauss; j ++) {
		//			int k;
		//			for(k = 0; k < n_int_bases[v]; k ++) {
		//				printf("%+.1e ",D[i][k][j]);
		//			} printf("\n");
		//		} printf("\n");
		//	} printf("\n");
		//	getchar();
		//}
	}

	//if(condition)
	//{
	//	for(i = 0; i < n_interpolations; i ++) {
	//		for(j = 0; j < n_gauss; j ++) {
	//			int k;
	//			for(k = 0; k < n_int_bases[interpolation_variable[i]]; k ++) {
	//				printf("%+.1e ",face->Q[i][k][j]);
	//			} printf("\n");
	//		} printf("\n");
	//	} printf("\n");
	//	getchar();
	//}

	return FACE_SUCCESS;
}

//////////////////////////////////////////////////////////////////
