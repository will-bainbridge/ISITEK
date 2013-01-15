#ifndef EXPRESSION_H
#define EXPRESSION_H

typedef struct s_EXPRESSION * EXPRESSION;

EXPRESSION expression_generate(char *string);
int expression_number_of_substitutes(EXPRESSION expression);
int expression_number_of_recursions(EXPRESSION expression);
void expression_evaluate(int n, double *value, EXPRESSION expression, double **substitute, double **work);
void expression_print(EXPRESSION expression);
void expression_free(EXPRESSION expression);

#endif
