///////////////////////////////////////////////////////////////////////////////

typedef struct s_EXPRESSION * EXPRESSION;

EXPRESSION expression_generate(char *string);
void expression_evaluate(int n, double *value, EXPRESSION expression, double **substitute, double **work);
void expression_destroy(EXPRESSION expression);

///////////////////////////////////////////////////////////////////////////////
