#ifndef CONDITION_H
#define CONDITION_H

typedef struct s_CONDITION * CONDITION;

CONDITION condition_empty();
CONDITION condition_condition(char *name);
char * condition_name(CONDITION condition);
int condition_n_variables(CONDITION condition);
int condition_max_n_variables();
void condition_variable(CONDITION condition, int *variable);
void condition_differential(CONDITION condition, int *differential);
int condition_n_parameters(CONDITION condition);
void (* condition_value(CONDITION condition))(double *, double *, double *, double *);

#endif
