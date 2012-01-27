//////////////////////////////////////////////////////////////////

#include<stdio.h>
#include<stdlib.h>

#include "isitek.h"

void read_geometry(FILE *file, int *n_nodes, struct NODE **node, int *n_faces, struct FACE **face, int *n_elements, struct ELEMENT **element);
void process_geometry(int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element);

void boundaries_input(FILE *file, int n_faces, struct FACE *face, int *n_boundaries, struct BOUNDARY **boundary);
void terms_input(FILE *file, int *n_terms, struct TERM **term);

void update_face_integration(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face);
void update_element_integration(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element);

void update_element_unknowns(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element, int n_u_old, int *n_u, double *u_old, double **u);

void update_face_boundaries(int n_variables, int n_faces, struct FACE *face, int n_boundaries, struct BOUNDARY *boundary);

void update_element_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element);
void update_face_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face, int n_boundaries_old, struct BOUNDARY *boundary_old);

void initialise_values(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, double *initial, double *u);
void initialise_system(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, int n_u, SPARSE *system);

void calculate_system(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, int n_terms, struct TERM *term, int n_u, double *u, SPARSE system, double *residual);

void write_case(FILE *file, int n_variables, int *variable_order, int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element, int n_boundaries, struct BOUNDARY *boundary);
void read_case(FILE *file, int *n_variables, int **variable_order, int *n_nodes, struct NODE **node, int *n_faces, struct FACE **face, int *n_elements, struct ELEMENT **element, int *n_boundaries, struct BOUNDARY **boundary);

void generate_numbered_file_path(char *file_path, char *base_path, int number);

void read_data(FILE *file, int *n_u, double **u, int *number);
void write_data(FILE *file, int n_u, double *u, int number);

//////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	exit_if_false(argc == 2,"exactly two input arguments required");

	int i, n;

	// open the input file
	char *input_file_path = argv[1];
	FILE *input_file = fopen(input_file_path,"r");
	exit_if_false(input_file != NULL,"opening input file");

	// allocate the file paths
	char *case_file_path, *data_file_path, *data_numbered_file_path;
	exit_if_false(case_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating case file path");
	exit_if_false(data_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating data file path");
	exit_if_false(data_numbered_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating data numbered file path");

	// get the case file path from the input file
	exit_if_false(fetch_value(input_file,"case_file_path",'s',case_file_path) == FETCH_SUCCESS,"reading case_file_path from the input file");

	// mesh structures
	int n_nodes, n_faces, n_elements, n_boundaries_old = 0, n_boundaries, n_terms;
	struct NODE *node;
	struct FACE *face;
	struct ELEMENT *element;
	struct BOUNDARY *boundary_old = NULL, *boundary;
	struct TERM *term;
	
	// solution vectors
	int n_u_old = 0, n_u;
	double *u_old = NULL, *u;

	// get the number of variables from the input file
	int n_variables_old = 0, n_variables;
	exit_if_false(fetch_value(input_file,"number_of_variables",'i',&n_variables) == FETCH_SUCCESS,"reading number_of_variables from the input file");

	// allocate initial values
	double *initial = (double *)malloc(n_variables * sizeof(double));
	exit_if_false(initial != NULL,"allocating initial values");

	// read the orders
	int *variable_order_old = NULL, *variable_order = (int *)malloc(n_variables * sizeof(int));
	exit_if_false(variable_order != NULL,"allocating orders");
	exit_if_false(fetch_vector(input_file, "variable_order", 'i', n_variables, variable_order) == FETCH_SUCCESS,"reading variable_order from the input file");

	// read the variable names
	char **variable_name = allocate_character_matrix(NULL,n_variables,MAX_STRING_LENGTH);
	exit_if_false(variable_name != NULL,"allocating variable names");
	warn_if_false(fetch_vector(input_file,"variable_name",'s',n_variables,variable_name) == FETCH_SUCCESS,"reading variable_name from the input file");

	// iteration counters
	int outer_iteration = 0, inner_iteration;
	int n_outer_iterations, n_inner_iterations;
	exit_if_false(fetch_value(input_file,"number_of_inner_iterations",'i',&n_inner_iterations) == FETCH_SUCCESS,"reading number_of_inner_iterations from the input file");
	exit_if_false(fetch_value(input_file,"number_of_outer_iterations",'i',&n_outer_iterations) == FETCH_SUCCESS,"reading number_of_outer_iterations from the input file");

	// open the case file
	FILE *case_file = fopen(case_file_path,"r");

	// read the case file
	if(case_file != NULL)
	{
		// read the case file
		read_case(case_file, &n_variables_old, &variable_order_old, &n_nodes, &node, &n_faces, &face, &n_elements, &element, &n_boundaries_old, &boundary_old);
		fclose(case_file);

		// numbers of uknowns
		n = 0; for(i = 0; i < n_variables; i ++) n += n_elements*ORDER_TO_N_BASIS(variable_order[i]);

		// read initial data file if it exists
		if(fetch_value(input_file,"initial_data_file_path",'s',data_file_path) == FETCH_SUCCESS)
		{
			FILE *initial_data_file = fopen(data_file_path,"r");
			exit_if_false(initial_data_file != NULL,"opening the initial data file");
			read_data(initial_data_file, &n_u_old, &u_old, &outer_iteration);
			fclose(initial_data_file);
			exit_if_false(n_u_old == n,"case and initial data files do not match");
		}
	}

	// read the geometry file
	else
	{
	 	// get geometry file path from input file
		char *geometry_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
		exit_if_false(geometry_file_path != NULL,"allocating geometry_file_path");
		exit_if_false(fetch_value(input_file,"geometry_file_path",'s',geometry_file_path) == FETCH_SUCCESS,"reading geometry_file_path from the input file");

		// open the geometry file
		FILE *geometry_file = fopen(geometry_file_path,"r");
		exit_if_false(geometry_file != NULL,"opening geometry file");
		free(geometry_file_path);

	 	// read and process geometry
		read_geometry(geometry_file, &n_nodes, &node, &n_faces, &face, &n_elements, &element);
		fclose(geometry_file);
		process_geometry(n_nodes, node, n_faces, face, n_elements, element);
	}

	// get the data file path
	exit_if_false(fetch_value(input_file,"data_file_path",'s',data_file_path) == FETCH_SUCCESS,"reading data_file_path from the input file");

	// read boundaries and equation terms from input_file
	boundaries_input(input_file, n_faces, face, &n_boundaries, &boundary);
	terms_input(input_file, &n_terms, &term);

	// update unknown indices and values
	update_element_unknowns(n_variables_old, n_variables, variable_order_old, variable_order, n_elements, element, n_u_old, &n_u, u_old, &u);

	// update face boundaries
	update_face_boundaries(n_variables, n_faces, face, n_boundaries, boundary);

	// update integration
	update_face_integration(n_variables_old, n_variables, variable_order_old, variable_order, n_faces, face);
	update_element_integration(n_variables_old, n_variables, variable_order_old, variable_order, n_elements, element);

	// update numerics
	update_face_numerics(n_variables_old, n_variables, variable_order_old, variable_order, n_faces, face, n_boundaries_old, boundary_old);
	update_element_numerics(n_variables_old, n_variables, variable_order_old, variable_order, n_elements, element);

	// initialise the values
	if(fetch_vector(input_file,"variable_initial_value",'d',n_variables,initial) == FETCH_SUCCESS)
		initialise_values(n_variables, variable_order, n_elements, element, initial, u);

	// initialise the system and residuals
	SPARSE system = NULL;
	initialise_system(n_variables, variable_order, n_elements, element, n_u, &system);
	double *residual = (double *)malloc(n_u * sizeof(double));
	exit_if_false(residual != NULL,"allocating the residuals");

	// iterate
	n_outer_iterations += outer_iteration;
	for(; outer_iteration < n_outer_iterations; outer_iteration ++)
	{
		for(inner_iteration = 0; inner_iteration < n_inner_iterations; inner_iteration ++)
		{
			// generate system
			calculate_system(n_variables, variable_order, n_elements, element, n_terms, term, n_u, u, system, residual);
			//

			// solve
		}

		generate_numbered_file_path(data_numbered_file_path, data_file_path, outer_iteration + 1);
		FILE *data_file = fopen(data_numbered_file_path,"w");
		exit_if_false(data_file != NULL,"opening data file");
		write_data(data_file, n_u, u, outer_iteration + 1);
		fclose(data_file);
	}

	// save to case file
	case_file = fopen(case_file_path,"w");
	exit_if_false(case_file != NULL,"opening case file");
	write_case(case_file, n_variables, variable_order, n_nodes, node, n_faces, face, n_elements, element, n_boundaries, boundary);
	fclose(case_file);

	// display

	// clean up
	fclose(input_file);

	free(case_file_path);
	free(data_file_path);
	free(data_numbered_file_path);

	destroy_nodes(n_nodes,node);
	destroy_faces(n_faces,face,n_variables);
	destroy_elements(n_elements,element,n_variables);
	destroy_boundaries(n_boundaries_old,boundary_old);
	destroy_boundaries(n_boundaries,boundary);
	destroy_terms(n_terms,term);

	free(initial);

	free(variable_order_old);
	free(variable_order);
	destroy_matrix((void *)variable_name);

	free(u_old);
	free(u);

	sparse_destroy(system);
	free(residual);

	return 0;
}

//////////////////////////////////////////////////////////////////
