//////////////////////////////////////////////////////////////////

#include<stdio.h>
#include<stdlib.h>

#include "isitek.h"

void read_geometry(FILE *file, int *n_nodes, struct NODE **node, int *n_faces, struct FACE **face, int *n_elements, struct ELEMENT **element);
void process_geometry(int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element);

void boundaries_input(FILE *file, int n_faces, struct FACE *face, int *n_boundaries, struct BOUNDARY **boundary);
void terms_input(FILE *file, int *n_terms, struct TERM **term);
int initial_input(FILE *file, int n_variables, EXPRESSION **initial);

void update_element_unknowns(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element, int n_u_old, int *n_u, double *u_old, double **u);
void update_face_boundaries(int n_variables, int n_faces, struct FACE *face, int n_boundaries, struct BOUNDARY *boundary);

int update_face_integration(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face);
int update_element_integration(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element);

int update_face_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_faces, struct FACE *face, int n_boundaries_old, struct BOUNDARY *boundary_old);
int update_element_numerics(int n_variables_old, int n_variables, int *variable_order_old, int *variable_order, int n_elements, struct ELEMENT *element);

void initialise_values(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, EXPRESSION *initial, double *u);
void initialise_system(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, int n_u, SPARSE *system);

void calculate_system(int n_variables, int *variable_order, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element, int n_terms, struct TERM *term, int n_u, double *u_old, double *u, SPARSE system, double *residual);
void slope_limit(int n_variables, int *variable_order, int n_nodes, struct NODE *node, int n_elements, struct ELEMENT *element, int n_boundaries, struct BOUNDARY *boundary, double *u);
void calculate_maximum_changes_and_residuals(int n_variables, int *variable_order, int n_elements, struct ELEMENT *element, double *du, double *max_u, double *residual, double *max_residual);

void write_case(FILE *file, int n_variables, int *variable_order, int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element, int n_boundaries, struct BOUNDARY *boundary);
void read_case(FILE *file, int *n_variables, int **variable_order, int *n_nodes, struct NODE **node, int *n_faces, struct FACE **face, int *n_elements, struct ELEMENT **element, int *n_boundaries, struct BOUNDARY **boundary);

void generate_numbered_file_path(char *file_path, char *base_path, int number);

void read_data(FILE *file, int *n_u, double **u, int *number);
void write_data(FILE *file, int n_u, double *u, int number);

void write_display(FILE *file, int n_variables, char **variable_name, int *variable_order, int n_elements, struct ELEMENT *element, int n_u, double *u);

//////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	exit_if_false(argc == 2,"exactly one input argument required");

	timer_start();

	//-------------------------------------------------------------------//
	
	// counters
	int i, n;

	// file paths
	char *input_file_path, *geometry_file_path, *case_file_path, *data_file_path, *data_numbered_file_path, *display_file_path, *display_numbered_file_path;
	exit_if_false(geometry_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating geometry file path");
	exit_if_false(case_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating case file path");
	exit_if_false(data_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating data file path");
	exit_if_false(data_numbered_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating data numbered file path");
	exit_if_false(display_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating display file path");
	exit_if_false(display_numbered_file_path = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating display numbered file path");

	// print string
	char *print;
	exit_if_false(print = (char *)malloc(MAX_STRING_LENGTH * sizeof(char)),"allocating the print string");

	// mesh structures
	int n_nodes, n_faces, n_elements, n_boundaries_old = 0, n_boundaries, n_terms;
	struct NODE *node;
	struct FACE *face;
	struct ELEMENT *element;
	struct BOUNDARY *boundary_old = NULL, *boundary;
	struct TERM *term;
	EXPRESSION *initial = NULL;
	
	// solution vectors
	int n_u_old = 0, n_u;
	double *u_old = NULL, *u;

	// files
	FILE *input_file, *case_file, *data_file, *geometry_file, *display_file;

	//-------------------------------------------------------------------//

	// opening the input file
	print_info("opening the input file %s",argv[1]);
	input_file_path = argv[1];
	input_file = fopen(input_file_path,"r");
	exit_if_false(input_file != NULL,"opening %s",input_file_path);

	// reading the case file path
	print_info("reading the case file path");
	exit_if_false(fetch_value(input_file,"case_file_path",'s',case_file_path) == FETCH_SUCCESS,"reading case_file_path from %s",input_file_path);
	print_continue(case_file_path);

	// reading the number of variables, variable names and orders
	print_info("reading the variables");
	int n_variables_old = 0, n_variables;
	exit_if_false(fetch_value(input_file,"number_of_variables",'i',&n_variables) == FETCH_SUCCESS,"reading number_of_variables from %s",input_file_path);
	int *variable_order_old = NULL, *variable_order = (int *)malloc(n_variables * sizeof(int));
	exit_if_false(variable_order != NULL,"allocating orders");
	exit_if_false(fetch_vector(input_file, "variable_order", 'i', n_variables, variable_order) == FETCH_SUCCESS,"reading variable_order from %s",input_file_path);
	char **variable_name = allocate_character_matrix(NULL,n_variables,MAX_STRING_LENGTH);
	exit_if_false(variable_name != NULL,"allocating variable names");
	warn_if_false(fetch_vector(input_file,"variable_name",'s',n_variables,variable_name) == FETCH_SUCCESS,"reading variable_name from %s",input_file_path);
	for(i = 0; i < n_variables; i ++) print_continue("%s order %i",variable_name[i],variable_order[i]);

	// reading the number of inner and outer iterations to perform
	print_info("reading the numbers of iterations");
	int outer_iteration = 0, inner_iteration;
	int n_outer_iterations, n_inner_iterations, data_n_outer_iterations, display_n_outer_iterations;
	exit_if_false(fetch_value(input_file,"number_of_inner_iterations",'i',&n_inner_iterations) == FETCH_SUCCESS,"reading number_of_inner_iterations from %s",input_file_path);
	exit_if_false(fetch_value(input_file,"number_of_outer_iterations",'i',&n_outer_iterations) == FETCH_SUCCESS,"reading number_of_outer_iterations from %s",input_file_path);
	print_continue("%i outer of %i inner iterations to be done",n_outer_iterations,n_inner_iterations);

	// read existing case and data
	case_file = fopen(case_file_path,"r");
	if(case_file != NULL)
	{
		print_info("reading existing case file %s",case_file_path);
		read_case(case_file, &n_variables_old, &variable_order_old, &n_nodes, &node, &n_faces, &face, &n_elements, &element, &n_boundaries_old, &boundary_old);
		fclose(case_file);
		n = 0; for(i = 0; i < n_variables; i ++) n += n_elements*ORDER_TO_N_BASIS(variable_order[i]);

		if(fetch_value(input_file,"initial_data_file_path",'s',data_file_path) == FETCH_SUCCESS)
		{
			print_info("reading existing data file %s",data_file_path);
			exit_if_false(data_file = fopen(data_file_path,"r"),"opening %s",data_file_path);
			read_data(data_file, &n_u_old, &u_old, &outer_iteration);
			fclose(data_file);
			exit_if_false(n_u_old == n,"case and initial data does not match");
		}
	}

	// construct new case
	else
	{
		print_info("reading the geometry file path");
		exit_if_false(fetch_value(input_file,"geometry_file_path",'s',geometry_file_path) == FETCH_SUCCESS,"reading geometry_file_path from %s",input_file_path);
		print_continue(geometry_file_path);

		print_info("reading the geometry file %s",geometry_file_path);
		exit_if_false(geometry_file = fopen(geometry_file_path,"r"),"opening %s",geometry_file_path);
		read_geometry(geometry_file, &n_nodes, &node, &n_faces, &face, &n_elements, &element);
		fclose(geometry_file);
		print_continue("%i nodes, %i faces and %i elements",n_nodes,n_faces,n_elements);

		print_info("generating additional connectivity and geometry");
		process_geometry(n_nodes, node, n_faces, face, n_elements, element);
	}

	// read the data file path and output frequency
	print_info("reading the data file path and output frequency");
	exit_if_false(fetch_value(input_file,"data_file_path",'s',data_file_path) == FETCH_SUCCESS,"reading data_file_path from %s",input_file_path);
	if(fetch_value(input_file,"data_number_of_outer_iterations",'i',&data_n_outer_iterations) != FETCH_SUCCESS)
		data_n_outer_iterations = n_outer_iterations + outer_iteration;
	print_continue("data to be written to %s every %i outer iterations",data_file_path,data_n_outer_iterations);

	// read boundaries
	print_info("reading boundaries");
	boundaries_input(input_file, n_faces, face, &n_boundaries, &boundary);
	print_continue("%i boundaries",n_boundaries);

	// read terms
	print_info("reading PDE terms");
	terms_input(input_file, &n_terms, &term);
	print_continue("%i terms",n_terms);

	// update unknown indices and values
	print_info("updating the numbering of the degrees of freedom");
	update_element_unknowns(n_variables_old, n_variables, variable_order_old, variable_order, n_elements, element, n_u_old, &n_u, u_old, &u);
	print_continue("%i degrees of freedom",n_u);

	// update face boundaries
	print_info("updating the face boundary associations");
	update_face_boundaries(n_variables, n_faces, face, n_boundaries, boundary);

	// update integration
	print_info("updating integration");
	i = update_face_integration(n_variables_old, n_variables, variable_order_old, variable_order, n_faces, face);
	if(i) print_continue("updated %i face quadratures",i);
	i = update_element_integration(n_variables_old, n_variables, variable_order_old, variable_order, n_elements, element);
	if(i) print_continue("updated %i element quadratures",i);

	// update numerics
	print_info("updating numerics");
	i = update_face_numerics(n_variables_old, n_variables, variable_order_old, variable_order, n_faces, face, n_boundaries_old, boundary_old);
	if(i) print_continue("updated %i face interpolations",i);
	i = update_element_numerics(n_variables_old, n_variables, variable_order_old, variable_order, n_elements, element);
	if(i) print_continue("updated %i element interpolations",i);

	// write case file
	print_info("writing case file %s",case_file_path);
	exit_if_false(case_file = fopen(case_file_path,"w"),"opening %s",case_file_path);
	write_case(case_file, n_variables, variable_order, n_nodes, node, n_faces, face, n_elements, element, n_boundaries, boundary);
	fclose(case_file);

	// read the display file path and output frequency
	print_info("reading the display file path and output frequency");
	if(
			fetch_value(input_file,"display_file_path",'s',display_file_path) == FETCH_SUCCESS &&
			fetch_value(input_file,"display_number_of_outer_iterations",'i',&display_n_outer_iterations) == FETCH_SUCCESS
	  )
		print_continue("display to be written to %s every %i outer iterations",display_file_path,display_n_outer_iterations);
	else
	{
		display_n_outer_iterations = 0;
		warn_if_false(0,"display files will not be written");
	}

	// initialise
	if(initial_input(input_file, n_variables, &initial))
	{
		print_info("initialising the degrees of freedom");
		initialise_values(n_variables, variable_order, n_elements, element, initial, u);
	}

	//-------------------------------------------------------------------//
	
	// allocate and initialise the system
	print_info("allocating and initialising the system");
	SPARSE system = NULL;
	initialise_system(n_variables, variable_order, n_elements, element, n_u, &system);
	double *residual, *max_residual, *du, *max_du;
	exit_if_false(residual = (double *)malloc(n_u * sizeof(double)),"allocating the residuals");
	exit_if_false(max_residual = (double *)malloc(n_variables * sizeof(double)),"allocating the maximum residuals");
	exit_if_false(du = (double *)malloc(n_u * sizeof(double)),"allocating du");
	exit_if_false(max_du = (double *)malloc(n_variables * sizeof(double)),"allocating the maximum changes");
	exit_if_false(u_old = (double *)realloc(u_old, n_u * sizeof(double)),"re-allocating u_old");

	timer_reset();

	// iterate
	print_info("iterating");
	n_outer_iterations += outer_iteration;
	for(; outer_iteration < n_outer_iterations; outer_iteration ++)
	{
		print_output("iteration %i", outer_iteration);

		for(i = 0; i < n_u; i ++) u_old[i] = u[i];

		for(inner_iteration = 0; inner_iteration < n_inner_iterations; inner_iteration ++)
		{
			calculate_system(n_variables, variable_order, n_faces, face, n_elements, element, n_terms, term, n_u, u_old, u, system, residual);

			exit_if_false(sparse_solve(system, du, residual) == SPARSE_SUCCESS,"solving system");
			for(i = 0; i < n_u; i ++) u[i] -= du[i];

			calculate_maximum_changes_and_residuals(n_variables, variable_order, n_elements, element, du, max_du, residual, max_residual);
			for(i = 0; i < n_variables; i ++) sprintf(&print[26*i],"%.6e|%.6e ",max_du[i],max_residual[i]);
			print_continue("%s",print);
		}

		//slope_limit(n_variables, variable_order, n_nodes, node, n_elements, element, n_boundaries, boundary, u);

		if(data_n_outer_iterations && outer_iteration % data_n_outer_iterations == 0)
		{
			generate_numbered_file_path(data_numbered_file_path, data_file_path, outer_iteration);
			print_info("writing data to %s",data_numbered_file_path);
			exit_if_false(data_file = fopen(data_numbered_file_path,"w"),"opening %s",data_numbered_file_path);
			write_data(data_file, n_u, u, outer_iteration);
			fclose(data_file);
		}

		if(display_n_outer_iterations && outer_iteration % display_n_outer_iterations == 0)
		{
			generate_numbered_file_path(display_numbered_file_path, display_file_path, outer_iteration);
			print_info("writing display to %s",display_numbered_file_path);
			exit_if_false(display_file = fopen(display_numbered_file_path,"w"),"opening %s",display_numbered_file_path);
			write_display(display_file, n_variables, variable_name, variable_order, n_elements, element, n_u, u);
			fclose(display_file);
		}

		timer_print();
	}

	//-------------------------------------------------------------------//
	
	print_info("freeing all memory");

	fclose(input_file);

	free(geometry_file_path);
	free(case_file_path);
	free(data_file_path);
	free(data_numbered_file_path);
	free(display_file_path);
	free(display_numbered_file_path);

	free(print);

	destroy_nodes(n_nodes,node);
	destroy_faces(n_faces,face,n_variables);
	destroy_elements(n_elements,element,n_variables);
	destroy_boundaries(n_boundaries_old,boundary_old);
	destroy_boundaries(n_boundaries,boundary);
	destroy_terms(n_terms,term);
	destroy_initial(n_variables,initial);

	free(variable_order_old);
	free(variable_order);
	destroy_matrix((void *)variable_name);

	free(u_old);
	free(u);

	sparse_destroy(system);
	free(residual);
	free(max_residual);
	free(du);
	free(max_du);

	return 0;
}

//////////////////////////////////////////////////////////////////
