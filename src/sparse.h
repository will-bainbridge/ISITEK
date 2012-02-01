#ifndef SPARSE_H
#define SPARSE_H

typedef struct s_SPARSE * SPARSE;

SPARSE sparse_allocate(SPARSE sparse, int n_rows);
SPARSE sparse_allocate_rows(SPARSE sparse, int *n_non_zeros);
void sparse_destroy(SPARSE sparse);
void sparse_print(SPARSE sparse);
void sparse_set_row_indices(SPARSE sparse, int row, int *index);
void sparse_set_row_values(SPARSE sparse, int row, double *value);
void sparse_add_to_row_values(SPARSE sparse, int row, double *value);
void sparse_set_zero(SPARSE sparse);
int sparse_solve_umfpack(SPARSE sparse, double *x, double *b);

#define SPARSE_SUCCESS 1
#define SPARSE_SOLVE_ERROR -1
#define SPARSE_MEMORY_ERROR -1

#endif

