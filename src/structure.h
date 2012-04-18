#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "expression.h"

struct NODE
{
	double x[2];
};

struct FACE
{
	int n_nodes;
	struct NODE **node;

	int n_borders;
	struct ELEMENT **border;

	int *n_boundaries;
	struct BOUNDARY ***boundary;

	double normal[2];
	double centre[2];
	double size;

	double   **X; // n_points * 2
	double    *W; // n_points
	double ****Q; // n_variables * n_differentials * ( n_borders * n_basis + n_boundaries ) * n_points
};

struct ELEMENT
{
	int n_faces;
	struct FACE **face;

	double centre[2];
	double size;

	int **unknown;

	double  **X; // n_points * 2
	double   *W; // n_points
	double ***P; // n_differentials * n_basis * n_points
	double ***Q; // n_faces * n_basis * n_points
	double ***I; // n_variables * n_points * n_basis
};

struct BOUNDARY
{
	int n_faces;
	struct FACE **face;

	int variable;
	int condition[2];
	double value;
};

struct TERM
{
	int equation;
	char type;
	double implicit;
	int n_variables;
	int *variable;
	int *differential;
	char *method;
	EXPRESSION *weight;
	EXPRESSION residual;
	EXPRESSION *jacobian;
};

#endif
