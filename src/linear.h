#ifndef LINEAR_H
#define LINEAR_H

// blas routines
void dcopy_(int const * n, double const * X, int const * incx, double * Y, int const * incy);
void dscal_(int const * n, double const * alpha, double * X, int const * incx);
double ddot_(int const * n, double const * X, int const * incx, double const * Y, int const * incy);
void daxpy_(int const * n, double const * alpha, double const * X, int const * incx, double * Y, int const * incy);
void dgemv_(char const * trans, int const * m, int const * n, double const * alpha, double const * A, int * const lda, double const * X, int const * incx, double const * beta, double * Y, int const * incy);
void dgemm_(char const * transa, char const * transb, int const * m, int const * n, int const * k, double const * alpha, double const * A, int const * lda, double const * B, int const * ldb, double const * beta, double * C, int const * ldc);

// lapack routines
void dgesv_(int const * n, int const * nrhs, double * A, int const * lda, int * ipiv, double * B, int const * ldb, int * info);

#endif
