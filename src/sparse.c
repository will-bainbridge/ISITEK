//////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "sparse.h"

//////////////////////////////////////////////////////////////////

struct s_SPARSE
{
	int n_rows; // number of rows
	int *row; // row indices or pointers
	int *column; // column indices
	double *value; // values

	int n_sub_matrices; // number of sub matrices
	int **sub_matrix_shape; // sub matrix numbers of rows and columns
	int **sub_matrix_index; // sub matrix indices
};

//////////////////////////////////////////////////////////////////

void sort(int number, int length, int **index, int *value);
void sort_sift_down(int number, int start, int end, int **index, int *value);
int sort_less_than(int number, int a, int b, int **index);
#define SORT_SWAP(a,b) do{ int t = a; a = b; b = t; } while(0)

//////////////////////////////////////////////////////////////////

/*int main()
{
	int anr = 2, ar[2] = {0,2};
	int anc = 2, ac[2] = {1,3};
	double ax[4] = {1,1,1,1};

	int bnr = 3, br[3] = {1,2,3};
	int bnc = 2, bc[2] = {0,1};
	double bx[6] = {2,2,2,2,2,2};

	SPARSE s = sparse_new(NULL);

	sparse_insert_sub_matrix(s,anr,ar,anc,ac);
	sparse_insert_sub_matrix(s,bnr,br,bnc,bc);
	sparse_insert_sub_matrix(s,bnc,bc,bnr,br);
	sparse_order_sub_matrices(s);
	sparse_add_sub_matrix(s,0,ax);
	sparse_add_sub_matrix(s,2,bx);
	sparse_add_sub_matrix(s,1,bx);
	sparse_print(s);
	sparse_spy(s);
	sparse_free(s);

	return 0;
}*/

//////////////////////////////////////////////////////////////////

SPARSE sparse_allocate(SPARSE sparse)
{
	SPARSE new;

	new = (SPARSE)realloc(sparse, sizeof(struct s_SPARSE));
	if(new == NULL) return NULL;

	new->n_rows = 0;
	new->n_sub_matrices = 0;

	if(sparse == NULL)
	{
		new->row = NULL;
		new->column = NULL;
		new->value = NULL;
		new->sub_matrix_shape = NULL;
		new->sub_matrix_index = NULL;
	}

	return new;
}

//////////////////////////////////////////////////////////////////

int sparse_insert_sub_matrix(SPARSE sparse, int n_rows, int *row, int n_columns, int *column)
{
	int i, j, n_indices = 0;

	sparse->n_sub_matrices += 1;

	int **new = (int **)realloc(sparse->sub_matrix_shape, sparse->n_sub_matrices * sizeof(int *));
	if(new == NULL) return SPARSE_MEMORY_ERROR;
	if(sparse->sub_matrix_shape == NULL) new[0] = NULL;
	new[0] = (int *)realloc(new[0], 2 * sparse->n_sub_matrices * sizeof(int));
	if(new[0] == NULL) return SPARSE_MEMORY_ERROR;
	for(i = 1; i < sparse->n_sub_matrices; i++) new[i] = new[i-1] + 2;
	sparse->sub_matrix_shape = new;

	sparse->sub_matrix_shape[sparse->n_sub_matrices-1][0] = n_rows;
	sparse->sub_matrix_shape[sparse->n_sub_matrices-1][1] = n_columns;

	for(i = 0; i < sparse->n_sub_matrices; i ++) n_indices += sparse->sub_matrix_shape[i][0]*sparse->sub_matrix_shape[i][1];

	sparse->row = (int *)realloc(sparse->row, n_indices * sizeof(int));
	if(sparse->row == NULL) return SPARSE_MEMORY_ERROR;

	sparse->column = (int *)realloc(sparse->column, n_indices * sizeof(int));
	if(sparse->column == NULL) return SPARSE_MEMORY_ERROR;

	n_indices -= n_rows*n_columns;

	for(i = 0; i < n_rows; i ++)
	{
		for(j = 0; j < n_columns; j ++)
		{
			sparse->row[n_indices] = row[i];
			sparse->column[n_indices] = column[j];
			n_indices ++;
		}
	}

	return sparse->n_sub_matrices - 1;;
}

//////////////////////////////////////////////////////////////////

int sparse_order_sub_matrices(SPARSE sparse)
{
	int i, j, n_indices = 0;

	// total number of indices
	for(i = 0; i < sparse->n_sub_matrices; i ++) n_indices += sparse->sub_matrix_shape[i][0]*sparse->sub_matrix_shape[i][1];

	// sort the indices
	int *number = (int *)malloc(n_indices * sizeof(int));
	if(number == NULL) return SPARSE_MEMORY_ERROR;
	for(i = 0; i < n_indices; i ++) number[i] = i;
	int *index[2] = { sparse->row , sparse->column };
	sort(2, n_indices, index, number);

	// allocate the sub matrix indices
	int **new = sparse->sub_matrix_index;
	sparse->sub_matrix_index = (int **)realloc(sparse->sub_matrix_index, sparse->n_sub_matrices * sizeof(int *));
	if(sparse->sub_matrix_index == NULL) return SPARSE_MEMORY_ERROR;
	if(new == NULL) sparse->sub_matrix_index[0] = NULL;
	sparse->sub_matrix_index[0] = (int *)realloc(sparse->sub_matrix_index[0], n_indices * sizeof(int));
	if(sparse->sub_matrix_index[0] == NULL) return SPARSE_MEMORY_ERROR;
	for(i = 0; i < sparse->n_sub_matrices - 1; i ++) sparse->sub_matrix_index[i+1] = sparse->sub_matrix_index[i] + sparse->sub_matrix_shape[i][0]*sparse->sub_matrix_shape[i][1];

	// label and remove duplicate sorted indices
	i = 0;
	j = 0;
	while(j < n_indices)
	{
		i += (sparse->row[i] != sparse->row[j]) || (sparse->column[i] != sparse->column[j]);
		sparse->sub_matrix_index[0][j] = i;
		sparse->row[i] = sparse->row[j];
		sparse->column[i] = sparse->column[j];
		j ++;
	}

	// number of non zeros
	int nnz = i + 1;

	// convert row indices to pointers
	// WARNING // this bit operates in-place and doesn't work if there are empty rows
	i = 0;
	j = 0;
	while(j < nnz - 1)
	{
		j ++;
		if(sparse->row[j-1] != sparse->row[j])
		{
			i ++;
			sparse->row[i] = j;
		}
	}
	sparse->row[++i] = nnz;

	// number of rows
	sparse->n_rows = i;

	// reallocate
	sparse->row = (int *)realloc(sparse->row, (sparse->n_rows + 1) * sizeof(int));
	if(sparse->row == NULL) return SPARSE_MEMORY_ERROR;
	sparse->column = (int *)realloc(sparse->column, sparse->row[sparse->n_rows] * sizeof(int));
	if(sparse->column == NULL) return SPARSE_MEMORY_ERROR;
	sparse->value = (double *)realloc(sparse->value, sparse->row[sparse->n_rows] * sizeof(double));
	if(sparse->value == NULL) return SPARSE_MEMORY_ERROR;

	// set to zero
	sparse_set_zero(sparse);

	// sort the numbering
	sort(1, n_indices, &number, sparse->sub_matrix_index[0]);
	
	// clean up
	free(number);

	/*// debug
	for(i = 0; i < sparse->n_rows + 1; i ++) printf("%i ",sparse->row[i]);
	printf("\n");
	for(i = 0; i < sparse->row[sparse->n_rows]; i ++) printf("%i ",sparse->column[i]);
	printf("\n\n");
	for(i = 0; i < sparse->n_sub_matrices; i ++)
	{
		for(j = 0; j < sparse->sub_matrix_shape[i][0]*sparse->sub_matrix_shape[i][1]; j ++)
		{
			printf("%i ",sparse->sub_matrix_index[i][j]);
		}
		printf("\n");
	}
	printf("\n");*/

	return SPARSE_SUCCESS;
}

//////////////////////////////////////////////////////////////////

void sparse_set_sub_matrix(SPARSE sparse, int index, double **value)
{
	int i, j, k = 0;
	for(i = 0; i < sparse->sub_matrix_shape[index][0]; i ++)
	{
		for(j = 0; j < sparse->sub_matrix_shape[index][1]; j ++)
		{
			sparse->value[sparse->sub_matrix_index[index][k++]] = value[i][j];
		}
	}
}

//////////////////////////////////////////////////////////////////

void sparse_add_sub_matrix(SPARSE sparse, int index, double **value)
{
	int i, j, k = 0;
	for(i = 0; i < sparse->sub_matrix_shape[index][0]; i ++)
	{
		for(j = 0; j < sparse->sub_matrix_shape[index][1]; j ++)
		{
			sparse->value[sparse->sub_matrix_index[index][k++]] += value[i][j];
		}
	}
}

//////////////////////////////////////////////////////////////////

void sparse_destroy(SPARSE sparse)
{
	free(sparse->row);
	free(sparse->column);
	free(sparse->value);
	if(sparse->sub_matrix_shape) free(sparse->sub_matrix_shape[0]);
	free(sparse->sub_matrix_shape);
	if(sparse->sub_matrix_index) free(sparse->sub_matrix_index[0]);
	free(sparse->sub_matrix_index);
	free(sparse);
}

//////////////////////////////////////////////////////////////////

void sparse_print(SPARSE sparse)
{
	int i, j;
	for(i = 0; i < sparse->n_rows; i ++)
	{
		for(j = sparse->row[i]; j < sparse->row[i+1]; j ++)
		{
			printf("%5i %5i %+.10e\n",i,sparse->column[j],sparse->value[j]);
		}
	}
}

//////////////////////////////////////////////////////////////////

int sparse_spy(SPARSE sparse, int height, int width)
{
	int i, j;

	int n_levels = 7;
	char level[][8] = { "\033[0m" , "\033[45m" , "\033[44m" , "\033[46m" , "\033[42m" , "\033[43m" , "\033[41m" };

	int **count = (int **)malloc(height * sizeof(int *));
	if(count == NULL) return SPARSE_MEMORY_ERROR;
	count[0] = (int *)malloc(height * width * sizeof(int));
	if(count[0] == NULL) return SPARSE_MEMORY_ERROR;
	for(i = 1; i < height; i ++) count[i] = count[i-1] + width;
	for(i = 0; i < height*width; i ++) count[0][i] = 0;

	int size[2] = { sparse->n_rows , 0 };
	for(i = 0; i < sparse->row[sparse->n_rows]; i ++) size[1] = size[1] > sparse->column[i] ? size[1] : sparse->column[i];

	for(i = 0; i < sparse->n_rows; i ++)
		for(j = sparse->row[i]; j < sparse->row[i+1]; j ++)
			count[i*height/size[0]][sparse->column[j]*width/size[1]] ++;

	int max_count = 0;
	for(i = 0; i < height*width; i ++) max_count = max_count > count[0][i] ? max_count : count[0][i];

	for(i = 0; i < height; i ++)
	{
		for(j = 0; j < width; j ++)
		{
			printf("%s  ",level[n_levels*count[i][j]/(max_count+1)]);
		}
		printf("%s\n",level[0]);
	}

	free(count[0]);
	free(count);

	return SPARSE_SUCCESS;
}

//////////////////////////////////////////////////////////////////

void sparse_set_zero(SPARSE sparse)
{
	int i;
	for(i = 0; i < sparse->row[sparse->n_rows]; i ++) sparse->value[i] = 0.0;
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

	status = umfpack_di_symbolic(sparse->n_rows, sparse->n_rows, sparse->row, sparse->column, sparse->value, &symbolic, control, info);
	if(status < 0) return SPARSE_SOLVE_ERROR;

	status = umfpack_di_numeric(sparse->row, sparse->column, sparse->value, symbolic, &numeric, control, info);
	if(status < 0) return SPARSE_SOLVE_ERROR;

	umfpack_di_free_symbolic(&symbolic);

	status = umfpack_di_solve(UMFPACK_At, sparse->row, sparse->column, sparse->value, x, b, numeric, control, info);
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
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n_rows), sparse->value, sparse->row, sparse->column, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if(error != 0) return SPARSE_SOLVE_ERROR;

	// numerical factorisation
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n_rows), sparse->value, sparse->row, sparse->column, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if(error != 0) return SPARSE_SOLVE_ERROR;

	// back substitution and iterative refinement
	phase = 33;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n_rows), sparse->value, sparse->row, sparse->column, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	if(error != 0) return SPARSE_SOLVE_ERROR;

	// termination and release of memory
	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(sparse->n_rows), &ddum, sparse->row, sparse->column, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

	return SPARSE_SUCCESS;
}

//////////////////////////////////////////////////////////////////

#else
#error solver definition not found or recognised
#endif

//////////////////////////////////////////////////////////////////

void sort(int number, int length, int **index, int *value)
{
	int i;
	int start, end;

	for(start = (length-2)/2; start >= 0; start --)
	{
		sort_sift_down(number,start,length,index,value);
	}

	for(end = length-1; end > 0; end --)
	{
		for(i = 0; i < number; i ++) SORT_SWAP(index[i][0],index[i][end]);
		SORT_SWAP(value[0],value[end]);
		sort_sift_down(number,0,end,index,value);
	}
}

//////////////////////////////////////////////////////////////////

void sort_sift_down(int number, int start, int end, int **index, int *value)
{
	int i, root = start, child;

	while(2*root + 1 < end)
	{
		child = 2*root + 1;

		child += (child + 1 < end) && sort_less_than(number,child,child+1,index);

		if(sort_less_than(number,root,child,index))
		{
			for(i = 0; i < number; i ++) SORT_SWAP(index[i][child],index[i][root]);
			SORT_SWAP(value[child],value[root]);
			root = child;
		}

		else return;
	}
}

//////////////////////////////////////////////////////////////////

int sort_less_than(int number, int a, int b, int **index)
{
	int i, r = 0;

	for(i = 0; i < number; i ++)
	{
		if(index[i][a] == index[i][b]) continue;
		else { r = index[i][a] < index[i][b]; break; }
	}

	return r;
}

//////////////////////////////////////////////////////////////////
