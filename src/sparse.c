//////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "sparse.h"
//#include "umfpack.h"

//////////////////////////////////////////////////////////////////

struct s_SPARSE
{
	int n; //number of rows
	int *row; //row pointers (length n+1)
	int *index; //column indices (length nnz)
	double *value; //values (length nnz)
};

//////////////////////////////////////////////////////////////////

SPARSE sparse_allocate(SPARSE sparse, int n_rows)
{
	SPARSE new;

	new = (SPARSE)realloc(sparse,sizeof(struct s_SPARSE));
	if(new == NULL) return NULL;

	if(sparse == NULL) { new->row = new->index = NULL; new->value = NULL; }

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

	sparse->value = (double *)realloc(sparse->value,n*sizeof(double));
	if(sparse->value == NULL) return NULL;

	return sparse;
}

//////////////////////////////////////////////////////////////////

void sparse_destroy(SPARSE sparse)
{
	free(sparse->row);
	free(sparse->index);
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
			printf("%5i %5i %+15.10e\n",i,sparse->index[j],sparse->value[j]);
		}
	}
}

//////////////////////////////////////////////////////////////////

void sparse_set_row(SPARSE sparse, int row, int *index, double *value)
{
	int i, r, n;

	r = sparse->row[row];
	n = sparse->row[row+1] - r;

	for(i = 0; i < n; i ++)
	{
		sparse->index[r+i] = index[i];
		sparse->value[r+i] = value[i];
	}
}

//////////////////////////////////////////////////////////////////

/*int sparse_solve_umfpack(SPARSE sparse, double *x, double *b)
{
	void *symbolic, *numeric;

	umfpack_di_symbolic(sparse->n, sparse->n, sparse->row, sparse->index, sparse->value, &symbolic, NULL, NULL);

	umfpack_di_numeric(sparse->row, sparse->index, sparse->value, symbolic, &numeric, NULL, NULL);
	umfpack_di_free_symbolic(&symbolic);

	umfpack_di_solve(UMFPACK_At, sparse->row, sparse->index, sparse->value, x, b, numeric, NULL, NULL);
	umfpack_di_free_numeric(&numeric);

	return SPARSE_SUCCESS;
}*/

//////////////////////////////////////////////////////////////////
