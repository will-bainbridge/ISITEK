//////////////////////////////////////////////////////////////////

#include "isitek.h"

#include "constants.h"
#include "linear.h"

void basis(int n, double *phi, double **x, double *origin, double size, int index, int *differential);
void transformation_matrix(int order, double **T, double **R);

//////////////////////////////////////////////////////////////////

int update_face_integration(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face)
{
	int f, g, i, v;

	int max_variable_order = 0, max_variable_order_old = 0;
	for(v = 0; v < n_variables; v ++) max_variable_order = MAX(max_variable_order,variable_order[v]);
	for(v = 0; v < n_variables_old; v ++) max_variable_order_old = MAX(max_variable_order_old,variable_order_old[v]);

	if(max_variable_order == max_variable_order_old) return 0;

	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order);

	double dx[2], size;

	for(f = 0; f < n_faces; f ++)
	{
		exit_if_false(face[f].X = allocate_face_x(&face[f], ORDER_TO_N_GAUSS(max_variable_order)),"allocating face X");
		exit_if_false(face[f].W = allocate_face_w(&face[f], ORDER_TO_N_GAUSS(max_variable_order)),"allocating face W");

		for(i = 0; i < 2; i ++) dx[i] = face[f].node[1]->x[i] - face[f].node[0]->x[i];

		size = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		for(g = 0; g < n_gauss; g ++)
		{
			for(i = 0; i < 2; i ++)
			{
				face[f].X[i][g] = 
					0.5 * (1.0 - gauss_locations[n_gauss-1][g]) * face[f].node[0]->x[i] + 
					0.5 * (1.0 + gauss_locations[n_gauss-1][g]) * face[f].node[1]->x[i];
			}

			face[f].W[g] = 0.5 * size * gauss_weights[n_gauss-1][g];
		}
	}

	return n_faces;
}

//////////////////////////////////////////////////////////////////

int update_element_integration(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element)
{
	int e, f, h, i, j, k, o, v;

	int max_variable_order = 0, max_variable_order_old = 0;
	for(v = 0; v < n_variables; v ++) max_variable_order = MAX(max_variable_order,variable_order[v]);
	for(v = 0; v < n_variables_old; v ++) max_variable_order_old = MAX(max_variable_order_old,variable_order_old[v]);

	if(max_variable_order == max_variable_order_old) return 0;

	int n_hammer = ORDER_TO_N_HAMMER(max_variable_order);
	double dx[2][2], size;

	for(e = 0; e < n_elements; e ++)
	{
		exit_if_false(element[e].X = allocate_element_x(&element[e],n_hammer*(element[e].n_faces - 2)),"allocating element X");
		exit_if_false(element[e].W = allocate_element_w(&element[e],n_hammer*(element[e].n_faces - 2)),"allocating element W");

		f = 0;

		for(i = 0; i < element[e].n_faces - 2; i ++)
		{
			// triangulate
			f ++;
			while(element[e].face[f]->node[0] == element[e].face[0]->node[0] || element[e].face[f]->node[1] == element[e].face[0]->node[0]) f ++;
			o = element[e].face[f]->border[0] == &element[e];
			for(j = 0; j < 2; j ++)
				for(k = 0; k < 2; k ++)
					dx[j][k] = element[e].face[f]->node[j == o]->x[k] - element[e].face[0]->node[0]->x[k];

			// hammer locations and weights
			size = dx[0][0]*dx[1][1] - dx[0][1]*dx[1][0];
			for(h = 0; h < n_hammer; h ++)
			{
				for(j = 0; j < 2; j ++)
				{
					element[e].X[j][i*n_hammer+h] = element[e].face[0]->node[0]->x[j] +
						hammer_locations[n_hammer-1][0][h]*dx[0][j] +
						hammer_locations[n_hammer-1][1][h]*dx[1][j];
				}
				element[e].W[i*n_hammer+h] = hammer_weights[n_hammer-1][h] * size;
			}
		}
	}

	return n_elements;
}

//////////////////////////////////////////////////////////////////

void basis(int n, double *phi, double **x, double *origin, double size, int index, int *differential)
{
	int i, j;

	int zero = 0;
	for(i = 0; i < 2; i ++) zero += taylor_powers[index][i] < differential[i];
	if(zero)
	{
		for(i = 0; i < n; i ++) phi[i] = 0.0;
		return;
	}

	int power = 0;
	for(i = 0; i < 2; i ++) power += taylor_powers[index][i];

	double constant = taylor_coefficients[index] / pow(size,power);
	for(i = 0; i < 2; i ++) constant *= factorial[taylor_powers[index][i]] / factorial[taylor_powers[index][i] - differential[i]];

	for(i = 0; i < n; i ++)
	{
		phi[i] = constant;
		for(j = 0; j < 2; j ++) phi[i] *= pow( x[j][i] - origin[j] , taylor_powers[index][j] - differential[j] );
	}
}

//////////////////////////////////////////////////////////////////

int update_element_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element)
{
	int e, i, j, v;

	// old and new maximum variable orders
	int max_variable_order = 0, max_variable_order_old = 0;
	for(v = 0; v < n_variables; v ++) max_variable_order = MAX(max_variable_order,variable_order[v]);
	for(v = 0; v < n_variables_old; v ++) max_variable_order_old = MAX(max_variable_order_old,variable_order_old[v]);

	// what needs updating
	int max_update = max_variable_order != max_variable_order_old, any_update = n_variables_old < n_variables;
	for(v = 0; v < MIN(n_variables_old,n_variables); v ++) any_update = any_update || (variable_order_old[v] != variable_order[v]);
	if(!any_update) return 0;

	// numbers of basis functions
	int *n_basis, max_n_basis = ORDER_TO_N_BASIS(max_variable_order);
	exit_if_false(n_basis = (int *)malloc(n_variables * sizeof(int)),"allocating n_basis");
        for(v = 0; v < n_variables; v ++) n_basis[v] = ORDER_TO_N_BASIS(variable_order[v]);

	// numbers of points
	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order), n_points;

	// no differential
	int no_differential[2] = {0,0};

	// working matrices
	int lds = max_n_basis, ldm = max_n_basis, ldd = max_n_basis, lda = max_n_basis;
        int sizem = max_n_basis*max_n_basis;
	double **S = allocate_double_matrix(NULL,(MAX_ELEMENT_N_FACES-2)*n_hammer,lds);
	double **M = allocate_double_matrix(NULL,max_n_basis,ldm);
        double **D = allocate_double_matrix(NULL,max_n_basis,ldd);
	double **A = allocate_double_matrix(NULL,max_n_basis,lda);
	double **X = allocate_double_matrix(NULL,2,MAX_ELEMENT_N_FACES);
	exit_if_false(S != NULL && M != NULL && A != NULL && X != NULL,"allocating working matrices");

	// lapack and blas
	int info, *pivot = (int *)malloc((max_n_basis + 2) * sizeof(int));
        exit_if_false(pivot != NULL,"allocating pivot");
        char trans[2] = "NT";
        int int_1 = 1;
        double dbl_0 = 0.0, dbl_1 = 1.0;

	for(e = 0; e < n_elements; e ++)
	{
		n_points = n_hammer*(element[e].n_faces - 2);

		if(max_update)
		{
			// interior matrices
			exit_if_false(element[e].P = allocate_element_p(&element[e],max_n_basis,n_points),"allocating element P");
			for(i = 0; i < max_n_basis; i ++)
				for(j = 0; j < max_n_basis; j ++)
					basis(n_points,element[e].P[i][j],element[e].X,element[e].centre,element[e].size,j,taylor_powers[i]);

			// face matrices
			exit_if_false(element[e].Q = allocate_element_q(&element[e],max_n_basis,n_gauss),"allocating element Q");
			for(i = 0; i < element[e].n_faces; i ++)
				for(j = 0; j < max_n_basis; j ++)
					basis(n_gauss,element[e].Q[i][j],element[e].face[i]->X,element[e].centre,element[e].size,j,no_differential);

			// corner matrix
			exit_if_false(element[e].V = allocate_element_v(&element[e],max_n_basis),"allocating element V");
			for(i = 0; i < element[e].n_faces; i ++)
				for(j = 0; j < 2; j ++)
					X[j][i] = element[e].face[i]->node[element[e].face[i]->border[0] != &element[e]]->x[j];
			for(i = 0; i < max_n_basis; i ++)
				basis(element[e].n_faces,element[e].V[i],X,element[e].centre,element[e].size,i,no_differential);
		}

		// mass matrix
                for(i = 0; i < max_n_basis; i ++) dcopy_(&n_points,element[e].P[powers_taylor[0][0]][i],&int_1,&S[0][i],&lds);
                for(i = 0; i < n_points; i ++) dscal_(&max_n_basis,&element[e].W[i],S[i],&int_1);
                dgemm_(&trans[0],&trans[0],&max_n_basis,&max_n_basis,&n_points,&dbl_1,S[0],&lds,element[e].P[powers_taylor[0][0]][0],&n_points,&dbl_0,M[0],&ldm);

		// initialise matrices
		exit_if_false(element[e].I = allocate_element_i(&element[e],n_variables,n_basis,n_points),"allocating element I");
                for(v = 0; v < n_variables; v ++)
                {
			if(!max_update && n_variables_old > v) if(variable_order_old[v] == variable_order[v]) continue;
			for(i = 0; i < n_basis[v]; i ++) dcopy_(&n_points,&S[0][i],&lds,&element[e].I[v][0][i],&n_basis[v]);
			dcopy_(&sizem,M[0],&int_1,A[0],&int_1);
			dgesv_(&n_basis[v],&n_points,A[0],&lda,pivot,element[e].I[v][0],&n_basis[v],&info);
		}

		// limiting matrices
		if(max_n_basis > 1)
		{
			// diffusion matrix
			for(i = 0; i < max_n_basis; i ++) dcopy_(&n_points,element[e].P[powers_taylor[1][0]][i],&int_1,&S[0][i],&lds);
			for(i = 0; i < n_points; i ++) dscal_(&max_n_basis,&element[e].W[i],S[i],&int_1);
			dgemm_(&trans[0],&trans[0],&max_n_basis,&max_n_basis,&n_points,&dbl_1,S[0],&lds,element[e].P[powers_taylor[1][0]][0],&n_points,&dbl_0,D[0],&ldd);
			for(i = 0; i < max_n_basis; i ++) dcopy_(&n_points,element[e].P[powers_taylor[0][1]][i],&int_1,&S[0][i],&lds);
			for(i = 0; i < n_points; i ++) dscal_(&max_n_basis,&element[e].W[i],S[i],&int_1);
			dgemm_(&trans[0],&trans[0],&max_n_basis,&max_n_basis,&n_points,&dbl_1,S[0],&lds,element[e].P[powers_taylor[0][1]][0],&n_points,&dbl_1,D[0],&ldd);
		}

		// limiting matrices
		exit_if_false(element[e].L = allocate_element_l(&element[e],n_variables,n_basis),"allocating element L");
		for(v = 0; v < n_variables; v ++)
		{
			if(n_variables_old > v) if(variable_order_old[v] == variable_order[v]) continue;
			if(variable_order[v] == 1) continue;
			dcopy_(&sizem,M[0],&int_1,A[0],&int_1);
			for(i = 0; i < n_basis[v]; i ++) dcopy_(&n_basis[v],D[i],&int_1,element[e].L[v][i],&int_1);
			dgesv_(&n_basis[v],&n_basis[v],A[0],&lda,pivot,element[e].L[v][0],&n_basis[v],&info);
		}
	}

	free(n_basis);
	destroy_matrix((void *)S);
	destroy_matrix((void *)M);
	destroy_matrix((void *)D);
	destroy_matrix((void *)A);
	destroy_matrix((void *)X);
	free(pivot);

	return n_elements;
}

//////////////////////////////////////////////////////////////////

void transformation_matrix(int order, double **T, double **R)
{
	int i, j, k, n = ORDER_TO_N_BASIS(order), row[2], col[2];

	for(i = 0; i < n; i ++) for(j = 0; j < n; j ++) T[i][j] = 0.0;
	T[0][0] = 1.0;

	for(i = 1; i < order; i ++)
	{
		for(j = 0; j < ORDER_TO_N_BASIS(i); j ++)
		{
			if(taylor_powers[j][0] + taylor_powers[j][1] == i - 1)
			{
				for(k = 0; k < ORDER_TO_N_BASIS(i); k ++)
				{
					row[0] = powers_taylor[taylor_powers[j][0]+1][taylor_powers[j][1]];
					row[1] = powers_taylor[taylor_powers[j][0]][taylor_powers[j][1]+1];

					col[0] = powers_taylor[taylor_powers[k][0]+1][taylor_powers[k][1]];
					col[1] = powers_taylor[taylor_powers[k][0]][taylor_powers[k][1]+1];

					T[row[0]][col[0]] += R[0][0]*T[j][k];
					T[row[0]][col[1]] += R[0][1]*T[j][k];

					if(taylor_powers[j][0]) continue;

					T[row[1]][col[0]] += R[1][0]*T[j][k];
					T[row[1]][col[1]] += R[1][1]*T[j][k];
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////

int update_face_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face, int n_boundaries_old, struct BOUNDARY *boundary_old)
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
				for(j = 0; j < 2; j ++)
					if(face[f].boundary[v][i]->condition[j] != boundary_old[face_boundary_old[f][v][i]].condition[j])
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
					basis(n_gauss,&A[i][n_adj_bases],y,face[f].centre,face[f].size,face_taylor[i],bnd[b]->condition);

				for(i = 0; i < n_gauss; i ++)
					for(j = 0; j < n_int_terms; j ++)
						B[j][i+b*n_gauss+n_adj_bases] = B[i+b*n_gauss+n_adj_bases][j] = (i+b*n_gauss+n_adj_bases) == j;
			}

			/*if(n_bnd)
			{
				printf("\nface %i variable %i\n\n",f,v);

				for(i = 0; i < n_gauss; i ++) { printf("(%lf %lf) ",face[f].X[0][i],face[f].X[1][i]); } printf("\n\n");

				for(i = 0; i < n_int_bases; i ++) { printf("(%i,%i) ",taylor_powers[face_taylor[i]][0],taylor_powers[face_taylor[i]][1]); } printf("\n\n");

				for(i = 0; i < n_int_bases; i ++) { for(j = 0; j < n_int_bases; j ++) { printf("%+.1e ",A[j][i]); } printf("\n"); }
				for(i = 0; i < n_int_bases*9-1; i ++) { printf("-"); } printf("\n");
				for(i = 0; i < n_int_bases; i ++) { for(j = 0; j < n_int_terms; j ++) { printf("%+.1e ",B[j][i]); } printf("\n"); } printf("\n");
				
				getchar();
			}*/

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

			/*if(n_bnd)
			{
				for(j = 0; j < n_gauss; j ++) {
					for(i = 0; i < n_int_terms; i ++) {
						printf("%+.1e ",D[0][i][j]);
					} printf("\n");
				} printf("\n");
				getchar();
			}*/

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
}

//////////////////////////////////////////////////////////////////
