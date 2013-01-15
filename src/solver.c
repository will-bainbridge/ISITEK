////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "numerics.h"
#include "solver.h"

#define MAX_STRING_LENGTH 128
#define MAX_N_VARIABLES 5

////////////////////////////////////////////////////////////////////////////////

#include "INS.h"

////////////////////////////////////////////////////////////////////////////////

void solver_residual(double *u, double **r) { residual(u,r); }
void solver_jacobian(double *u, double ***j) { jacobian(u,j); }

static int variable_max_order = 0;
static int variable_n_bases[] = { 0 , 0 , 0 };
static int variable_max_n_bases = 0;
static int variable_sum_n_bases = 0;

static int n_gauss = 0;
static int n_hammer = 0;

////////////////////////////////////////////////////////////////////////////////

void solver_initialise()
{
	int i;
	for(i = 0; i < n_variables; i ++)
	{
		variable_max_order = variable_max_order > variable_order[i] ? variable_max_order : variable_order[i];
		variable_n_bases[i] = numerics_n_bases(variable_order[i]);
		variable_max_n_bases = variable_max_n_bases > variable_n_bases[i] ? variable_max_n_bases : variable_n_bases[i];
		variable_sum_n_bases += variable_n_bases[i];
	}
	n_gauss = numerics_n_gauss(variable_max_order);
	n_hammer = numerics_n_hammer(variable_max_order);
}

int solver_n_variables() { return n_variables; }
void solver_variable_order(int *order) { int i; for(i = 0; i < n_variables; i ++) order[i] = variable_order[i]; }
int solver_variable_max_order() { return variable_max_order; }
void solver_variable_n_bases(int *n_bases) { int i; for(i = 0; i < n_variables; i ++) n_bases[i] = variable_n_bases[i]; }
int solver_variable_max_n_bases() { return variable_max_n_bases; }
int solver_variable_sum_n_bases() { return variable_sum_n_bases; }

int solver_n_gauss() { return n_gauss; }
int solver_n_hammer() { return n_hammer; }
int solver_n_interpolations() { return n_interpolations; }
void solver_interpolation_variable(int *variable) { int i; for(i = 0; i < n_interpolations; i ++) variable[i] = interpolation_variable[i]; }
void solver_interpolation_differential(int *differential) { int i; for(i = 0; i < n_interpolations; i ++) differential[i] = interpolation_differential[i]; }
void solver_interpolation_method(char *method) { int i; for(i = 0; i < n_interpolations; i ++) method[i] = interpolation_method[i]; }

int solver_n_constants() { return n_constants; }
int solver_constant_set_value(char *name, double value)
{
	int i;
	for(i = 0; i < n_constants; i ++)
	{
		if(strcmp(name,constant_name[i]) == 0)
		{
			constant_value[i] = value;
			constant_set[i] = 1;
			break;
		}
	}
	return i < n_constants;
}
int solver_constants_set()
{
	int i, set = 1;
	for(i = 0; i < n_constants; i ++) set *= constant_set[i];
	return set;
}
int solver_constants_print(char *string)
{
	int i;
	for(i = 0; i < n_constants; i ++) if(sprintf(&string[strlen(string)],"%s=%g;",constant_name[i],constant_value[i]) < 0) break;
	return i == n_constants;
}

////////////////////////////////////////////////////////////////////////////////
