////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <string.h>
#include "condition.h"

#define MAX_STRING_LENGTH 128
#define MAX_N_CONSTRAINTS 5
#define MAX_N_VARIABLES 5
#define MAX_VARIABLE_N_CONSTRAINTS 2

////////////////////////////////////////////////////////////////////////////////

struct s_CONDITION
{
	char name[MAX_STRING_LENGTH];
	int n_constraints;
	int variable[MAX_N_CONSTRAINTS];
	int differential[MAX_N_CONSTRAINTS];
	int n_parameters;
	void (*value)(double *, double *, double *, double *);
	// void (*jacobian[])(double *, double *, double *, double *) // ???
	
	int variable_n_constraints[MAX_N_VARIABLES];
	int * variable_differential[MAX_N_VARIABLES], variable_differential_[MAX_N_VARIABLES * MAX_VARIABLE_N_CONSTRAINTS]; // should probably allocate these
};

////////////////////////////////////////////////////////////////////////////////

#include "INS_bc.h"

////////////////////////////////////////////////////////////////////////////////

static int max_n_constraints = 0;

////////////////////////////////////////////////////////////////////////////////

void condition_start()
{
	int i, j;

	for(i = 0; i < n_conditions; i ++)
		max_n_constraints = max_n_constraints > condition[i].n_constraints ? max_n_constraints : condition[i].n_constraints;

	for(i = 0; i < n_conditions; i ++)
		for(j = 0; j < MAX_N_VARIABLES; j ++)
		{
			condition[i].variable_n_constraints[j] = 0;
			condition[i].variable_differential[j] = condition[i].variable_differential_ + j*MAX_VARIABLE_N_CONSTRAINTS; // faff instead of dynamic allocation
		}

	for(i = 0; i < n_conditions; i ++)
		for(j = 0; j < condition[i].n_constraints; j ++)
			condition[i].variable_differential[condition[i].variable[j]][condition[i].variable_n_constraints[condition[i].variable[j]]++] =
				condition[i].differential[j];
}

//----------------------------------------------------------------------------//

CONDITION condition_empty() { return empty; }

//----------------------------------------------------------------------------//

CONDITION condition_condition(char *name)
{
	int i;
	for(i = 0; i < n_conditions; i ++) if(strcmp(name,condition[i].name) == 0) break;
	return i < n_conditions ? &condition[i] : NULL;
}

const char * condition_name(CONDITION condition)
{
	return condition->name;
}

//----------------------------------------------------------------------------//

int condition_n_constraints(CONDITION condition) { return condition->n_constraints; }
int condition_max_n_constraints() { return max_n_constraints; }
const int * condition_variable(CONDITION condition) { return condition->variable; }
const int * condition_differential(CONDITION condition) { return condition->differential; }
const int * condition_variable_n_constraints(CONDITION condition) { return condition->variable_n_constraints; }
const int * const * condition_variable_differential(CONDITION condition) { return condition->variable_differential; }

//----------------------------------------------------------------------------//

int condition_n_parameters(CONDITION condition) { return condition->n_parameters; }

//----------------------------------------------------------------------------//

//void (* condition_value(CONDITION condition))(double *, double *, double *, double *) { return condition->value; }

////////////////////////////////////////////////////////////////////////////////
