////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "numerics.h"
#include "solver.h"

////////////////////////////////////////////////////////////////////////////////

#include "INS.h"

////////////////////////////////////////////////////////////////////////////////

static int variable_max_order = 0;
static int *variable_n_bases;
static int variable_max_n_bases = 0;
static int variable_sum_n_bases = 0;

static int n_gauss = 0;
static int n_hammer = 0;

////////////////////////////////////////////////////////////////////////////////

int solver_start()
{
	int i;

	variable_n_bases = (int *)malloc(n_variables * sizeof(int));
	if(variable_n_bases == NULL) return SOLVER_MEMORY_ERROR;

	for(i = 0; i < n_variables; i ++)
	{
		variable_max_order = variable_max_order > variable_order[i] ? variable_max_order : variable_order[i];
		variable_n_bases[i] = numerics_n_bases(variable_order[i]);
		variable_max_n_bases = variable_max_n_bases > variable_n_bases[i] ? variable_max_n_bases : variable_n_bases[i];
		variable_sum_n_bases += variable_n_bases[i];
	}

	n_gauss = numerics_n_gauss(variable_max_order);
	n_hammer = numerics_n_hammer(variable_max_order);

	return SOLVER_SUCCESS;
}

void solver_end()
{
	free(variable_n_bases);
}

//----------------------------------------------------------------------------//

int solver_n_variables() { return n_variables; }
const int * solver_variable_order() { return variable_order; }
int solver_variable_max_order() { return variable_max_order; }
const int * solver_variable_n_bases() { return variable_n_bases; }
int solver_variable_max_n_bases() { return variable_max_n_bases; }
int solver_variable_sum_n_bases() { return variable_sum_n_bases; }

//----------------------------------------------------------------------------//

int solver_n_gauss() { return n_gauss; }
int solver_n_hammer() { return n_hammer; }

//----------------------------------------------------------------------------//

int solver_n_interpolations() { return n_interpolations; }
const int * solver_interpolation_variable() { return interpolation_variable; }
const int * solver_interpolation_differential() { return interpolation_differential; }
const char * solver_interpolation_method() { return interpolation_method; }

//----------------------------------------------------------------------------//

int solver_n_constants()
{
	return n_constants;
}

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

//----------------------------------------------------------------------------//

//void solver_residual(double *u, double **r) { residual(u,r); }
//void solver_jacobian(double *u, double ***j) { jacobian(u,j); }

////////////////////////////////////////////////////////////////////////////////
