#ifndef LINEAR_H
#define LINEAR_H

// blas routines
void dcopy_(int *n, double *X, int *incx, double *Y, int *incy);
void dscal_(int *n, double *alpha, double *X, int *incx);
void daxpy_(int *n, double *alpha, double *X, int *incx, double *Y, int *incy);
void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *lda, double *X, int *incx, double *beta, double *Y, int *incy);
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta, double *C, int *ldc);

// lapack routines
void dgesv_(int *n, int *nrhs, double *A, int *lda, int *ipiv, double *B, int *ldb, int *info);

#endif
