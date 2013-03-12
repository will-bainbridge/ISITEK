////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "boundary.h"
#include "condition.h"
#include "memory.h"

////////////////////////////////////////////////////////////////////////////////

struct s_BOUNDARY
{
	int index;
	char *name;
	CONDITION condition;
	int n_parameters;
	double *parameter;
};

////////////////////////////////////////////////////////////////////////////////

BOUNDARY boundary_new(int index)
{
	BOUNDARY boundary = (BOUNDARY)malloc(sizeof(struct s_BOUNDARY));
	if(boundary == NULL) return NULL;

	boundary->index = index;

	boundary->name = NULL;

	boundary->condition = condition_empty();

	boundary->n_parameters = 0;
	boundary->parameter = NULL;

	return boundary;
}

////////////////////////////////////////////////////////////////////////////////

int boundary_index(BOUNDARY boundary)
{
	return boundary->index;
}

////////////////////////////////////////////////////////////////////////////////

int boundary_read_name(FILE *file, BOUNDARY boundary)
{
	boundary->name = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	if(boundary->name == NULL) return BOUNDARY_MEMORY_ERROR;

	if(fscanf(file,"%s",boundary->name) != 1) return BOUNDARY_READ_ERROR;

	return BOUNDARY_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int boundary_read_faces(FILE *file, FACE *face, BOUNDARY boundary)
{
	int i, j, n, r, f[2];

	if(fscanf(file,"%i",&n) != 1) return BOUNDARY_READ_ERROR;

	for(i = 0; i < n; i ++)
	{
		r = fscanf(file,"%i:%i",&f[0],&f[1]);

		if(r < 1) return BOUNDARY_READ_ERROR;

		for(j = f[0]; j < (r==1)*(f[0]+1) + (r==2)*f[1]; j ++)
		{
			if(face_set_boundary(face[j],boundary) != FACE_SUCCESS) return BOUNDARY_LOGIC_ERROR;
		}
	}

	return BOUNDARY_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int boundary_read_condition(FILE *file, BOUNDARY boundary)
{
	int i;

	char *input = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	if(input == NULL) return BOUNDARY_MEMORY_ERROR;

	if(fscanf(file,"%s",input) != 1) return BOUNDARY_READ_ERROR;

	boundary->condition = condition_condition(input);
	if(boundary->condition == NULL) return BOUNDARY_READ_ERROR;

	boundary->n_parameters = condition_n_parameters(boundary->condition);

	boundary->parameter = (double *)malloc(boundary->n_parameters * sizeof(double));
	if(boundary->parameter == NULL) return BOUNDARY_MEMORY_ERROR;

	for(i = 0; i < boundary->n_parameters; i ++)
	{
		if(fscanf(file,"%lf",&(boundary->parameter[i])) != 1) return BOUNDARY_READ_ERROR;
	}

	free(input);

	return BOUNDARY_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

CONDITION boundary_condition(BOUNDARY boundary)
{
	return boundary->condition;
}

////////////////////////////////////////////////////////////////////////////////

void boundary_print(BOUNDARY boundary)
{
	printf("boundary %i\n",boundary->index);

	int i;
	if(boundary->name) printf("    boundary->name\n        %s\n",boundary->name);

	printf("    boundary->condition->name\n        %s\n",condition_name(boundary->condition));

	if(boundary->parameter)
	{
		printf("    boundary->parameter\n       ");
		for(i = 0; i < boundary->n_parameters; i ++) printf(" %g",boundary->parameter[i]);
		printf("\n");
	}
}

////////////////////////////////////////////////////////////////////////////////

void boundary_free(BOUNDARY boundary)
{
	free(boundary->name);
	free(boundary->parameter);
	free(boundary);
}

////////////////////////////////////////////////////////////////////////////////
