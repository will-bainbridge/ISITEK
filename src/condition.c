////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <string.h>
#include "condition.h"

#define MAX_STRING_LENGTH 128
#define MAX_N_VARIABLES 5

////////////////////////////////////////////////////////////////////////////////

struct s_CONDITION
{
	char name[MAX_STRING_LENGTH];
	int n_variables;
	int variable[MAX_N_VARIABLES];
	int differential[MAX_N_VARIABLES];
	int n_parameters;
	void (*value)(double *, double *, double *, double *);
	// void (*jacobian[])(double *, double *, double *, double *) // ???
};

////////////////////////////////////////////////////////////////////////////////

#include "INS_bc.h"

////////////////////////////////////////////////////////////////////////////////

CONDITION condition_empty() { return empty; }
CONDITION condition_condition(char *name)
{
	int i;
	for(i = 0; i < n_conditions; i ++) if(strcmp(name,condition[i].name) == 0) break;
	return i < n_conditions ? &condition[i] : NULL;
}

char * condition_name(CONDITION condition) { return condition->name; }
int condition_n_variables(CONDITION condition) { return condition->n_variables; }
int condition_max_n_variables()
{
	int i, max_n_variables = 0;
	for(i = 0; i < n_conditions; i ++) max_n_variables = max_n_variables > condition[i].n_variables ? max_n_variables : condition[i].n_variables;
	return max_n_variables;
}
void condition_variable(CONDITION condition, int *variable) { int i; for(i = 0; i < condition->n_variables; i ++) variable[i] = condition->variable[i]; }
void condition_differential(CONDITION condition, int *differential) { int i; for(i = 0; i < condition->n_variables; i ++) differential[i] = condition->differential[i]; }
int condition_n_parameters(CONDITION condition) { return condition->n_parameters; }
void (* condition_value(CONDITION condition))(double *, double *, double *, double *) { return condition->value; }

////////////////////////////////////////////////////////////////////////////////
