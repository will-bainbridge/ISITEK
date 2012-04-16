///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "expression.h"

#define LARGER '>'
#define SMALLER '<'
#define PLUS '+'
#define MINUS '-'
#define MULTIPLY '*'
#define DIVIDE '/'
#define POWER '^'

#define LEFTBRACE '('
#define RIGHTBRACE ')'

#define EQUALITY '='
#define SEPERATOR ';'
#define SUBSTITUTE '$'

#define SUBSTITUTE_ZERO '`' //must have a larger ascii value than all the others
#define VALUE '#'
#define EMPTY ' '
#define END '\0'

#define IS_OPERATOR(x) ( \
		(x) == LARGER || \
		(x) == SMALLER || \
		(x) == PLUS || \
		(x) == MINUS || \
		(x) == MULTIPLY || \
		(x) == DIVIDE || \
		(x) == POWER )

#define IS_CONTROL(x) ( \
		(x) == LEFTBRACE || \
		(x) == RIGHTBRACE )

#define IS_VALUE(x) ( \
		(x) == '0' || (x) == '1' || (x) == '2' || (x) == '3' || (x) == '4' || \
		(x) == '5' || (x) == '6' || (x) == '7' || (x) == '8' || (x) == '9' || \
		(x) == '.' )

#define IS_VARIABLE(x) ( \
		! IS_OPERATOR(x) && \
		! IS_CONTROL(x) && \
		! IS_VALUE(x) && \
		(x) != EQUALITY && \
		(x) != SEPERATOR && \
		(x) != SUBSTITUTE )

#define MAX_ELEMENTS 128
#define MAX_STRING 1024

static int precedence[128] = {
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,3,2,0,2,0,3,0,0, // * + - /
	0,0,0,0,0,0,0,0,0,0,
	1,0,1,0,0,0,0,0,0,0, // < >
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,4,0,0,0,0,0, // ^
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0};

struct s_EXPRESSION
{
	int type;
	double value;
};

void expression_simplify(EXPRESSION expression);

///////////////////////////////////////////////////////////////////////////////

/*int main()
{
	char *string = (char *)malloc(MAX_STRING * sizeof(char));
	strcpy(string,"cp=1005.0;cv=718.0;gamma=cp/cv;R=cp-cv;p0=1.0e5;((2.0*gamma*R/(gamma-1.0)*$1)*((p0/$0)^(gamma/(gamma-1.0))-1.0))^0.5");

	printf("    STRING > %s\n",string);

	EXPRESSION expression = expression_generate(string);
	if(expression == NULL) { printf("error generating the expression\n"); exit(1); }

	printf("EXPRESSION > ");
	expression_print(expression);
	printf("\n");

	int i;
	int n = 3;
	int n_s = expression_number_of_substitutes(expression);
	int n_r = expression_number_of_recursions(expression);

	printf("    N-SUBS > %i\n",n_s);
	printf("    N-RECS > %i\n",n_r);

	double *value = (double *)malloc(n * sizeof(double));
	double **substitute;
	substitute = (double **)malloc(n_s * sizeof(double *));
	substitute[0] = (double *)malloc(n_s * n * sizeof(double));
	for(i = 0; i < n_s; i ++) substitute[i] = substitute[0] + n*i;
	double **work;
	work = (double **)malloc(n_r * sizeof(double *));
	work[0] = (double *)malloc(n_r * n * sizeof(double));
	for(i = 0; i < n_r; i ++) work[i] = work[0] + n*i;

	substitute[0][0] = 9.5e4; substitute[0][1] = 9.6e4; substitute[0][2] = 9.6e4;
	substitute[1][0] = 300.0; substitute[1][1] = 300.0; substitute[1][2] = 310.0;

	expression_evaluate(n, value, expression, substitute, work);
	
	printf("    VALUES > %lf\n",value[0]);
	for(i = 1; i < n; i ++) printf("           > %lf\n",value[i]);

	free(string);
	expression_destroy(expression);

	free(value);
	free(substitute[0]);
	free(substitute);
	free(work[0]);
	free(work);

	return 0;
}*/

///////////////////////////////////////////////////////////////////////////////

EXPRESSION expression_generate(char *original)
{
	int i, length, offset;
	char *string = (char *)malloc(MAX_STRING * sizeof(char));
	strcpy(string,original);

	// remove whitespace

	offset = 0;
	length = strlen(string);

	for(i = 0; i < length; i ++)
	{
		if(string[i + offset] == ' ')
		{
			offset ++;
			length --;
		}
		string[i] = string[i + offset];
	}
	string[i] = '\0';

	// convert multiple expressions into just one

	length = strlen(string);

	int n_name, n_expr;
	char *name = (char *)malloc(MAX_STRING * sizeof(char));
	char *expr = (char *)malloc(MAX_STRING * sizeof(char));
	char *temp = (char *)malloc(MAX_STRING * sizeof(char));
	if(name == NULL || expr == NULL || temp == NULL) return NULL;

	while(IS_VARIABLE(string[0]))
	{
		i = 0;

		do name[i] = string[i]; while(string[++i] != EQUALITY);
		n_name = i;
		name[n_name] = '\0';

		i ++;

		do expr[i - n_name - 1] = string[i]; while(string[++i] != SEPERATOR);
		n_expr = i - n_name - 1;
		expr[n_expr] = '\0';

		i ++;

		while(string[i] == ';') i ++;

		sprintf(temp,"%s",&string[i]);
		strcpy(string,temp);

		i = 0;

		do
		{
			if(strncmp(&string[i],name,n_name) == 0 && (!IS_VARIABLE(string[i+n_name]) || string[i+n_name] == '\0'))
			{
				string[i] = '\0';

				sprintf(temp,"%s%c%s%c%s",string,LEFTBRACE,expr,RIGHTBRACE,&string[i + n_name]);

				strcpy(string,temp);
			}

			i ++;

		} while(i < strlen(string));
	}

	free(name);
	free(expr);
	free(temp);

	// convert the infix string into lists of postfix/RPN operations

	length = strlen(string);

	int o = 0;
	int *operator = (int *)malloc(MAX_ELEMENTS * sizeof(int));
	if(operator == NULL) return NULL;

	int e = 0;
	EXPRESSION expression = (EXPRESSION)malloc(MAX_ELEMENTS * sizeof(struct s_EXPRESSION));
	if(expression == NULL) return NULL;

	int index;

	i = 0;
	do {
		if(string[i] == LEFTBRACE)
		{
			operator[o++] = LEFTBRACE;
			i ++;
		}
		else if(string[i] == RIGHTBRACE)
		{
			while(operator[o-1] != LEFTBRACE)
			{
				expression[e++].type = operator[--o];
				if(!o) return NULL;
			}
			i ++;
			o --;
		}
		else if(IS_OPERATOR(string[i]))
		{
			while(o > 0)
				if(precedence[(int)string[i]] <= precedence[operator[o-1]]) expression[e++].type = operator[--o];
				else break;
			operator[o++] = string[i++];
		}
		else if(IS_VALUE(string[i]))
		{
			sscanf(&string[i],"%lf",&expression[e].value);
			expression[e++].type = VALUE;

			while((i < length) &&
					((!IS_OPERATOR(string[i]) && !IS_CONTROL(string[i])) ||
					 ((string[i] == MINUS || string[i] == PLUS) &&
					  (string[i-(i>0)] == 'e' || string[i-(i>0)] == 'E')))) i ++;
		}
		else if(string[i] == SUBSTITUTE)
		{
			sscanf(&string[++i],"%i",&index);
			expression[e++].type = SUBSTITUTE_ZERO + index;

			while(IS_VALUE(string[i])) i ++;
		}
		else return NULL;

	} while(i < length);

	while(o) expression[e++].type = operator[--o];

	expression[e].type = END;

	// simplify the expression as much as possible
	
	expression_simplify(expression);

	// clean up and return
	
	free(string);
	free(operator);

	e = 0;
	while(expression[e].type != END) e ++;

	return (EXPRESSION)realloc(expression, (e + 1) * sizeof(struct s_EXPRESSION));
}

///////////////////////////////////////////////////////////////////////////////

int expression_number_of_substitutes(EXPRESSION expression)
{
	int e = 0, n = 0;

	while(expression[e].type != END)
	{
		if(expression[e].type >= SUBSTITUTE_ZERO)
			n = n > expression[e].type - SUBSTITUTE_ZERO ? n : expression[e].type - SUBSTITUTE_ZERO;

		e ++;
	}

	return n + 1;
}

///////////////////////////////////////////////////////////////////////////////

int expression_number_of_recursions(EXPRESSION expression)
{
	int e = 0, w = 0, n = 0;

	while(expression[e].type != END)
	{
		if(expression[e].type == VALUE || expression[e].type >= SUBSTITUTE_ZERO)
			w ++;
		else if(IS_OPERATOR(expression[e].type))
			w --;

		n = n > w ? n : w;

		e ++;
	}

	return n;
}

///////////////////////////////////////////////////////////////////////////////

void expression_evaluate(int n, double *value, EXPRESSION expression, double **substitute, double **work)
{
	int i, e = 0, w = 0;

	while(expression[e].type != END)
	{
		if(expression[e].type == VALUE)
		{
			for(i = 0; i < n; i ++) work[w][i] = expression[e].value;
			w ++;
		}
		else if(expression[e].type >= SUBSTITUTE_ZERO)
		{
			for(i = 0; i < n; i ++) work[w][i] = substitute[expression[e].type - SUBSTITUTE_ZERO][i];
			w ++;
		}
		else if(IS_OPERATOR(expression[e].type))
		{
			switch(expression[e].type)
			{
				case POWER   : for(i = 0; i < n; i ++) work[w-2][i] = pow( work[w-2][i] , work[w-1][i] ); break;
				case MULTIPLY: for(i = 0; i < n; i ++) work[w-2][i] =      work[w-2][i] * work[w-1][i]  ; break;
				case DIVIDE  : for(i = 0; i < n; i ++) work[w-2][i] =      work[w-2][i] / work[w-1][i]  ; break;
				case PLUS    : for(i = 0; i < n; i ++) work[w-2][i] =      work[w-2][i] + work[w-1][i]  ; break;
				case MINUS   : for(i = 0; i < n; i ++) work[w-2][i] =      work[w-2][i] - work[w-1][i]  ; break;
				case LARGER  : for(i = 0; i < n; i ++) work[w-2][i] =      work[w-2][i] > work[w-1][i]  ; break;
				case SMALLER : for(i = 0; i < n; i ++) work[w-2][i] =      work[w-2][i] < work[w-1][i]  ; break;
			}
			w --;
		}
		e ++;
	}

	for(i = 0; i < n; i ++) value[i] = work[0][i];
}

///////////////////////////////////////////////////////////////////////////////

void expression_simplify(EXPRESSION expression)
{
	if(expression[0].type == END || expression[1].type == END) return;

	int i = 0, j, n;
	int type;
	double value;

	// fudge to make a non-allocated 2x1 matrix act like an allocated one
	double work[2], *p_work[2];
	p_work[0] = &work[0];
	p_work[1] = &work[1];

	// loop until a pass is completed with no simplifications
	do
	{
		// start at the begining haveing done no simplifications
		i = n = 0;

		while(1)
		{
			// find the next value-value-operator pattern
			while(!(expression[i].type == VALUE && expression[i+1].type == VALUE && IS_OPERATOR(expression[i+2].type)) && expression[i+2].type != END) i ++;

			// reached the end
			if(expression[i+2].type == END) break;

			// evaluate the sub-expression
			type = expression[i+3].type;
			expression[i+3].type = END;
			expression_evaluate(1, &value, &expression[i], NULL, p_work);
			expression[i+3].type = type;

			// replace the sub-expression with the value
			expression[i].value = value;
			j = i + 3;
			while(expression[j-1].type != END) { expression[j-2] = expression[j]; j ++; }

			// increment the number of simplifications
			n ++;
		}

	} while(n);
}

///////////////////////////////////////////////////////////////////////////////

void expression_print(EXPRESSION expression)
{
	int e = 0;
	while(expression[e].type != END)
	{
		if(IS_OPERATOR(expression[e].type) || IS_CONTROL(expression[e].type))
		{
			printf("%c",expression[e].type);
		}
		else if(expression[e].type == VALUE)
		{
			printf("{%g}",expression[e].value);
		}
		else if(expression[e].type >= SUBSTITUTE_ZERO)
		{
			printf("{$%i}",expression[e].type - SUBSTITUTE_ZERO);
		}
		e ++;
	}
}

///////////////////////////////////////////////////////////////////////////////

void expression_destroy(EXPRESSION expression)
{
	free(expression);
}

///////////////////////////////////////////////////////////////////////////////
