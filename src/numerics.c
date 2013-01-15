//////////////////////////////////////////////////////////////////

#include <math.h>
#include "numerics.h"

//////////////////////////////////////////////////////////////////

static double taylor_coefficient[55] = {1.0000000000000000e+00,1.0000000000000000e+00,1.0000000000000000e+00,5.0000000000000000e-01,1.0000000000000000e+00,5.0000000000000000e-01,1.6666666666666666e-01,4.9999999999999994e-01,4.9999999999999994e-01,1.6666666666666666e-01,4.1666666666666664e-02,1.6666666666666671e-01,2.5000000000000000e-01,1.6666666666666671e-01,4.1666666666666664e-02,8.3333333333333332e-03,4.1666666666666650e-02,8.3333333333333315e-02,8.3333333333333315e-02,4.1666666666666650e-02,8.3333333333333332e-03,1.3888888888888889e-03,8.3333333333333367e-03,2.0833333333333332e-02,2.7777777777777797e-02,2.0833333333333332e-02,8.3333333333333367e-03,1.3888888888888889e-03,1.9841269841269841e-04,1.3888888888888892e-03,4.1666666666666692e-03,6.9444444444444501e-03,6.9444444444444501e-03,4.1666666666666692e-03,1.3888888888888892e-03,1.9841269841269841e-04,2.4801587301587302e-05,1.9841269841269847e-04,6.9444444444444469e-04,1.3888888888888894e-03,1.7361111111111132e-03,1.3888888888888907e-03,6.9444444444444469e-04,1.9841269841269847e-04,2.4801587301587302e-05,2.7557319223985893e-06,2.4801587301587271e-05,9.9206349206349111e-05,2.3148148148148144e-04,3.4722222222222191e-04,3.4722222222222224e-04,2.3148148148148144e-04,9.9206349206349111e-05,2.4801587301587271e-05,2.7557319223985893e-06};

static int taylor_power[55][2] = {{0,0},{1,0},{0,1},{2,0},{1,1},{0,2},{3,0},{2,1},{1,2},{0,3},{4,0},{3,1},{2,2},{1,3},{0,4},{5,0},{4,1},{3,2},{2,3},{1,4},{0,5},{6,0},{5,1},{4,2},{3,3},{2,4},{1,5},{0,6},{7,0},{6,1},{5,2},{4,3},{3,4},{2,5},{1,6},{0,7},{8,0},{7,1},{6,2},{5,3},{4,4},{3,5},{2,6},{1,7},{0,8},{9,0},{8,1},{7,2},{6,3},{5,4},{4,5},{3,6},{2,7},{1,8},{0,9}};

static int power_taylor[10][10] = {
	{ 0, 2, 5, 9,14,20,27,35,44,54},
	{ 1, 4, 8,13,19,26,34,43,53, 0},
	{ 3, 7,12,18,25,33,42,52, 0, 0},
	{ 6,11,17,24,32,41,51, 0, 0, 0},
	{10,16,23,31,40,50, 0, 0, 0, 0},
	{15,22,30,39,49, 0, 0, 0, 0, 0},
	{21,29,38,48, 0, 0, 0, 0, 0, 0},
	{28,37,47, 0, 0, 0, 0, 0, 0, 0},
	{36,46, 0, 0, 0, 0, 0, 0, 0, 0},
	{45, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

static double factorial[10] = {1.0000000000000000e+00,1.0000000000000000e+00,2.0000000000000000e+00,6.0000000000000000e+00,2.4000000000000000e+01,1.2000000000000000e+02,7.2000000000000000e+02,5.0400000000000000e+03,4.0320000000000000e+04,3.6288000000000000e+05};

//////////////////////////////////////////////////////////////////

double numerics_taylor_coefficient(int index)
{
	return taylor_coefficient[index];
}

int numerics_taylor_power(int index, int dimension)
{
	return taylor_power[index][dimension];
}

int numerics_power_taylor(int power_x, int power_y)
{
	return power_taylor[power_x][power_y];
}

double numerics_factorial(int x)
{
	return factorial[x];
}

//////////////////////////////////////////////////////////////////

int numerics_n_bases(int order)
{
	return order*(order+1)/2;
}

int numerics_n_gauss(int order)
{
	return order;
}

int numerics_n_hammer(int order)
{
	return 12;
}

//////////////////////////////////////////////////////////////////

void numerics_basis(int n, double *phi, double **x, double *origin, double size, int index, int *differential)
{
	int i, j;

	int zero = 0;
	for(i = 0; i < 2; i ++) zero += taylor_power[index][i] < differential[i];
	if(zero)
	{
		for(i = 0; i < n; i ++) phi[i] = 0.0;
		return;
	}

	int power = 0;
	for(i = 0; i < 2; i ++) power += taylor_power[index][i];

	double constant = taylor_coefficient[index] / pow(size,power);
	for(i = 0; i < 2; i ++) constant *= factorial[taylor_power[index][i]] / factorial[taylor_power[index][i] - differential[i]];

	for(i = 0; i < n; i ++)
	{
		phi[i] = constant;
		for(j = 0; j < 2; j ++) phi[i] *= pow( x[j][i] - origin[j] , taylor_power[index][j] - differential[j] );
	}
}

//////////////////////////////////////////////////////////////////

void numerics_transformation_matrix(int order, double **T, double **R)
{
	int i, j, k, n = numerics_n_bases(order), row[2], col[2];

	for(i = 0; i < n; i ++) for(j = 0; j < n; j ++) T[i][j] = 0.0;
	T[0][0] = 1.0;

	for(i = 1; i < order; i ++)
	{
		for(j = 0; j < numerics_n_bases(i); j ++)
		{
			if(taylor_power[j][0] + taylor_power[j][1] == i - 1)
			{
				for(k = 0; k < numerics_n_bases(i); k ++)
				{
					row[0] = power_taylor[taylor_power[j][0]+1][taylor_power[j][1]];
					row[1] = power_taylor[taylor_power[j][0]][taylor_power[j][1]+1];

					col[0] = power_taylor[taylor_power[k][0]+1][taylor_power[k][1]];
					col[1] = power_taylor[taylor_power[k][0]][taylor_power[k][1]+1];

					T[row[0]][col[0]] += R[0][0]*T[j][k];
					T[row[0]][col[1]] += R[0][1]*T[j][k];

					if(taylor_power[j][0]) continue;

					T[row[1]][col[0]] += R[1][0]*T[j][k];
					T[row[1]][col[1]] += R[1][1]*T[j][k];
				}
			}
		}
	}
}

//////////////////////////////////////////////////////////////////

