//////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sparse.h"

//////////////////////////////////////////////////////////////////

struct s_SPARSE
{
	int n; //number of rows
	int *row; //row pointers (length n+1)
	int *index; //column indices (length nnz)
	int *order; //column insert order (length nnz)
	double *value; //values (length nnz)
};

void heap_sort(int *index, int *value, int n);
void sift_down(int *index, int *value, int lower, int upper);

//////////////////////////////////////////////////////////////////

SPARSE sparse_allocate(SPARSE sparse, int n_rows)
{
	SPARSE new;

	new = (SPARSE)realloc(sparse,sizeof(struct s_SPARSE));
	if(new == NULL) return NULL;

	if(sparse == NULL) { new->row = new->index = new->order = NULL; new->value = NULL; }

	new->n = n_rows;

	new->row = (int *)realloc(new->row,(n_rows+1)*sizeof(int));
	if(new->row == NULL) return NULL;

	return new;
}

//////////////////////////////////////////////////////////////////

SPARSE sparse_allocate_rows(SPARSE sparse, int *n_non_zeros)
{
	int i;

	sparse->row[0] = 0;
	for(i = 0; i < sparse->n; i ++) sparse->row[i+1] = sparse->row[i] + n_non_zeros[i];

	int n = sparse->row[sparse->n];

	sparse->index = (int *)realloc(sparse->index,n*sizeof(int));
	if(sparse->index == NULL) return NULL;

	sparse->order = (int *)realloc(sparse->order,n*sizeof(int));
	if(sparse->order == NULL) return NULL;

	sparse->value = (double *)realloc(sparse->value,n*sizeof(double));
	if(sparse->value == NULL) return NULL;

	return sparse;
}

//////////////////////////////////////////////////////////////////

void sparse_destroy(SPARSE sparse)
{
	free(sparse->row);
	free(sparse->index);
	free(sparse->order);
	free(sparse->value);
	free(sparse);
}

//////////////////////////////////////////////////////////////////

void sparse_print(SPARSE sparse)
{
	int i, j;
	for(i = 0; i < sparse->n; i ++)
	{
		for(j = sparse->row[i]; j < sparse->row[i+1]; j ++)
		{
			printf("%5i %5i %+.10e\n",i,sparse->index[j],sparse->value[j]);
		}
	}
}

//////////////////////////////////////////////////////////////////

void sparse_set_row_indices(SPARSE sparse, int row, int *index)
{
	int i, r = sparse->row[row], n = sparse->row[row+1] - r;

	for(i = 0; i < n; i ++) sparse->order[r+i] = i;

	heap_sort(index, &sparse->order[r], n);

	for(i = 0; i < n; i ++) sparse->index[r+i] = index[i];

	for(i = 0; i < n; i ++) index[sparse->order[r+i]] = sparse->index[r+i];
}

//////////////////////////////////////////////////////////////////

void sparse_set_row_values(SPARSE sparse, int row, double *value)
{
	int i, r = sparse->row[row], n = sparse->row[row+1] - r;

	for(i = 0; i < n; i ++) sparse->value[r+i] = value[sparse->order[r+i]];
}

//////////////////////////////////////////////////////////////////

void sparse_add_to_row_values(SPARSE sparse, int row, double *value)
{
	int i, r = sparse->row[row], n = sparse->row[row+1] - r;

	for(i = 0; i < n; i ++) sparse->value[r+i] += value[sparse->order[r+i]];
}

//////////////////////////////////////////////////////////////////

void sparse_set_zero(SPARSE sparse)
{
	int i;
	for(i = 0; i < sparse->row[sparse->n]; i ++) sparse->value[i] = 0.0;
}

//////////////////////////////////////////////////////////////////

#if defined(SOLVE_UMFPACK)

#include "umfpack.h"

int sparse_solve(SPARSE sparse, double *x, double *b)
{
	double info[UMFPACK_INFO], control[UMFPACK_CONTROL];
	void *symbolic, *numeric;
	int status;

	umfpack_di_defaults(control);

	status = umfpack_di_symbolic(sparse->n, sparse->n, sparse->row, sparse->index, sparse->value, &symbolic, control, info);
	if(status < 0) return SPARSE_SOLVE_ERROR;

	status = umfpack_di_numeric(sparse->row, sparse->index, sparse->value, symbolic, &numeric, control, info);
	if(status < 0) return SPARSE_SOLVE_ERROR;

	umfpack_di_free_symbolic(&symbolic);

	status = umfpack_di_solve(UMFPACK_At, sparse->row, sparse->index, sparse->value, x, b, numeric, control, info);
	if(status < 0) return SPARSE_SOLVE_ERROR;

	umfpack_di_free_numeric(&numeric);

	return SPARSE_SUCCESS;
}

//////////////////////////////////////////////////////////////////

#elif defined(SOLVE_PARDISO)

#include "mkl_pardiso.h"
#include "mkl_types.h"

int sparse_solve(SPARSE sparse, double *x, double *b)
{
	// auxilliary variables
	int i, idum;
	double ddum;

	// control parameters
	int maxfct, mnum, mtype, phase, nrhs, iparm[64], msglvl, error;
	maxfct = 1; // max 1 numerical factorisations
	mnum = 1; // sngle matrix to factorise
	mtype = 11; // real unsymmetric matrix
	nrhs = 1; // single right hand side
	for(i = 0; i < 64; i++) iparm[i] = 0; // initialise integer parameters
	iparm[0] = 1; // no defaults
	iparm[1] = 2; // metis reordering
	iparm[2] = 1; // 1 processor
	iparm[3] = 0; // no iterative-direct algorithm
	iparm[4] = 0; // no user fill-in reducing permutation
	iparm[5] = 0; // solution in x
	iparm[7] = 2; // max 2 iterative refinement steps
	iparm[9] = 13; // perturb the pivot elements with 1E-13
	iparm[10] = 1; // nonsymmetric permutation and scaling MPS
	iparm[12] = 1; // maximum weighted matching algorithm
	iparm[34] = 1; // c-style indexing from 0
	msglvl = 0; // print nothing
	error = 0; // initialise error flag

	// internal solver memory pointer
	void *pt[64];
	for(i = 0; i < 64; i++) pt[i] = 0;

	// reordering and symbolic factorisation
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n), sparse->value, sparse->row, sparse->index, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if(error != 0) return SPARSE_SOLVE_ERROR;

	// numerical factorisation
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n), sparse->value, sparse->row, sparse->index, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if(error != 0) return SPARSE_SOLVE_ERROR;

	// back substitution and iterative refinement
	phase = 33;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n), sparse->value, sparse->row, sparse->index, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if(error != 0) return SPARSE_SOLVE_ERROR;

	// termination and release of memory
	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n), &ddum, sparse->row, sparse->index, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	return SPARSE_SUCCESS;
}

//////////////////////////////////////////////////////////////////

#else
#error solver definition not found or recognised
#endif

//////////////////////////////////////////////////////////////////

void heap_sort(int *index, int *value, int n)
{
	int i;
	int itemp, vtemp;

	for(i = n/2; i >= 1; i--) sift_down(index - 1, value - 1, i, n);

	for(i = n; i >= 2; i--)
	{
		itemp = index[0];      vtemp = value[0];
		index[0] = index[i-1]; value[0] = value[i-1];
		index[i-1] = itemp;    value[i-1] = vtemp;

		sift_down(index - 1, value - 1, 1, i - 1);
	}
}

//////////////////////////////////////////////////////////////////

void sift_down(int *index, int *value, int lower, int upper)
{
	int i = lower, c = lower, lastindex = upper/2;
	int itemp, vtemp;

	itemp = index[i]; vtemp = value[i];

	while(c <= lastindex)
	{
		c = 2*i;
		if((c + 1 <= upper) && (index[c + 1] > index[c])) c++;

		if(itemp >= index[c]) break;

		index[i] = index[c]; value[i] = value[c];
		i = c;
	}

	index[i] = itemp; value[i] = vtemp;
}

//////////////////////////////////////////////////////////////////
