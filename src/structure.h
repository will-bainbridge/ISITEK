#ifndef STRUCTURE_H
#define STRUCTURE_H

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
	int *orient;

	double centre[2];
	double size;

	int **unknown;

	double   **X; // n_points * 2
	double    *W; // n_points
	double  ***P; // n_differentials * n_basis * n_points
	double  ***Q; // n_faces * n_basis * n_points
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
	int n_variables;
	int *variable;
	int *differential;
	int *power;
	char *method;
	double implicit;
	double constant;
};

#endif
