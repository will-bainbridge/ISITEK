#ifndef NUMERICS_H
#define NUMERICS_H

double numerics_taylor_coefficient(int index);
int numerics_taylor_power(int index, int dimension);
int numerics_power_taylor(int power_x, int power_y);
double numerics_factorial(int x);

int numerics_n_bases(int order);
int numerics_n_gauss(int order);
int numerics_n_hammer(int order);

void numerics_basis(int n, double *phi, double **x, double *origin, double size, int index, int *differential);
void numerics_transformation_matrix(int order, double **T, double **R);

#endif
