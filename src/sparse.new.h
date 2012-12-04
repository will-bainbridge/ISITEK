#ifndef SPARSE_H
#define SPARSE_H

typedef struct s_SPARSE * SPARSE;

SPARSE sparse_allocate(SPARSE sparse);
int sparse_insert_sub_matrix(SPARSE sparse, int n_rows, int *row, int n_columns, int *column);
int sparse_order_sub_matrices(SPARSE sparse);
void sparse_set_sub_matrix(SPARSE sparse, int index, double *value);
void sparse_add_sub_matrix(SPARSE sparse, int index, double *value);
void sparse_set_zero(SPARSE sparse);
void sparse_print(SPARSE sparse);
int sparse_spy(SPARSE sparse, int height, int width);
int sparse_solve(SPARSE sparse, double *x, double *b);
void sparse_destroy(SPARSE sparse);

#define SPARSE_SUCCESS 1
#define SPARSE_SOLVE_ERROR -1
#define SPARSE_MEMORY_ERROR -2

#endif

