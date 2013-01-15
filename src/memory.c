////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include "memory.h"

////////////////////////////////////////////////////////////////////////////////

int ** matrix_integer_new(int **old, int height, int width)
{
	int i;

	int **new = (int **)realloc(old, height * sizeof(int *));
	if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

	new[0] = (int *)realloc(new[0], height * width * sizeof(int));
	if(new[0] == NULL) return NULL;

	for(i = 1; i < height; i++) new[i] = new[i-1] + width;

	return new;
}

double ** matrix_double_new(double **old, int height, int width)
{
	int i;

	double **new = (double **)realloc(old, height * sizeof(double *));
	if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

	new[0] = (double *)realloc(new[0], height * width * sizeof(double));
	if(new[0] == NULL) return NULL;

	for(i = 1; i < height; i++) new[i] = new[i-1] + width;

	return new;
}

char ** matrix_character_new(char **old, int height, int width)
{
	int i;

	char **new = (char **)realloc(old, height * sizeof(char *));
	if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

	new[0] = (char *)realloc(new[0], height * width * sizeof(char));
	if(new[0] == NULL) return NULL;

	for(i = 1; i < height; i++) new[i] = new[i-1] + width;

	return new;
}

////////////////////////////////////////////////////////////////////////////////

void matrix_free(void **matrix)
{
	if(matrix)
	{
		free(matrix[0]);
	}
	free(matrix);
}

////////////////////////////////////////////////////////////////////////////////

int *** tensor_integer_new(int ***old, int height, int width, int depth)
{
        int i, j;

        int ***new = (int ***)realloc(old, height * sizeof(int **));
        if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

        new[0] = (int **)realloc(new[0], height * width * sizeof(int *));
        if(new[0] == NULL) return NULL;
	if(old == NULL) new[0][0] = NULL;

        new[0][0] = (int *)realloc(new[0][0], height * width * depth * sizeof(int));
        if(new[0][0] == NULL) return NULL;

        for(i = 1; i < height; i ++) new[i] = new[i-1] + width;
        for(i = 1; i < height; i ++) new[i][0] = new[i-1][0] + width*depth;
        for(i = 0; i < height; i ++) for(j = 1; j < width; j ++) new[i][j] = new[i][j-1] + depth;

        return new;
}

double *** tensor_double_new(double ***old, int height, int width, int depth)
{
        int i, j;

        double ***new = (double ***)realloc(old, height * sizeof(double **));
        if(new == NULL) return NULL;
	if(old == NULL) new[0] = NULL;

        new[0] = (double **)realloc(new[0], height * width * sizeof(double *));
        if(new[0] == NULL) return NULL;
	if(old == NULL) new[0][0] = NULL;

        new[0][0] = (double *)realloc(new[0][0], height * width * depth * sizeof(double));
        if(new[0][0] == NULL) return NULL;

        for(i = 1; i < height; i ++) new[i] = new[i-1] + width;
        for(i = 1; i < height; i ++) new[i][0] = new[i-1][0] + width*depth;
        for(i = 0; i < height; i ++) for(j = 1; j < width; j ++) new[i][j] = new[i][j-1] + depth;

        return new;
}

////////////////////////////////////////////////////////////////////////////////

void tensor_free(void ***tensor)
{
        if(tensor)
	{
		if(tensor[0])
		{
			free(tensor[0][0]);
		}
		free(tensor[0]);
	}
        free(tensor);
}

////////////////////////////////////////////////////////////////////////////////
