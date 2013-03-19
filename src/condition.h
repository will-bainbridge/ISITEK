#ifndef CONDITION_H
#define CONDITION_H

typedef struct s_CONDITION * CONDITION;

void condition_start();
CONDITION condition_empty();
CONDITION condition_condition(char *name);
const char * condition_name(CONDITION condition);
int condition_n_constraints(CONDITION condition);
int condition_max_n_constraints();
const int * condition_variable(CONDITION condition);
const int * condition_differential(CONDITION condition);
const int * condition_variable_n_constraints(CONDITION condition);
const int * const * condition_variable_differential(CONDITION condition);
int condition_n_parameters(CONDITION condition);
void (* condition_value(CONDITION condition))(double *, double *, double *, double *);

#endif
