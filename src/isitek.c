//////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "boundary.h"
#include "condition.h"
#include "element.h"
#include "expression.h"
#include "face.h"
#include "info.h"
#include "memory.h"
#include "node.h"
#include "solver.h"
#include "sparse.h"
#include "timer.h"
#include "utility.h"

//////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
	// counters
	int i, j;

	// structures
	int n_nodes, n_faces, n_elements, n_boundaries;
	NODE *node = NULL;
	FACE *face = NULL;
	ELEMENT *element = NULL;
	BOUNDARY *boundary = NULL;

	// initialisation
	EXPRESSION *initial;

	//-------------------------------------------------------------------//
	
	// initialise the solver
	exit_if_false(solver_start() == SOLVER_SUCCESS,"initialising the solver");
	condition_start();
	
	//-------------------------------------------------------------------//

	// opening the input file
	char *input_file_path = argv[1];
	FILE *input_file = fopen(input_file_path,"r");
	exit_if_false(input_file != NULL,"opening %s",input_file_path);

	// input data
	double value;
	char *string = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	exit_if_false(string != NULL,"allocating input string");

	// data positions
	fpos_t position[6] = {{-1},{-1},{-1},{-1},{-1},{-1}};
	while(fscanf(input_file,"%s",string) == 1)
	{
		i =
			1*(strcmp(string,"NODES") == 0) + 
			2*(strcmp(string,"EDGES") == 0 || strcmp(string,"FACES") == 0) + 
			3*(strcmp(string,"CELLS") == 0 || strcmp(string,"ELEMENTS") == 0) +
			4*(strcmp(string,"BOUNDARIES") == 0) + 
			5*(strcmp(string,"CONSTANTS") == 0) +
			6*(strcmp(string,"INITIAL") == 0);

		if(i) exit_if_false(fgetpos(input_file, &position[i-1]) == 0,"setting data position");
	}

	// read the nodes
	exit_if_false(fsetpos(input_file,&position[0]) == 0,"finding node data");
	exit_if_false(fscanf(input_file,"%i",&n_nodes) == 1,"reading the number of nodes");
	exit_if_false(node = (NODE *)malloc(n_nodes * sizeof(NODE)),"allocating nodes");
	for(i = 0; i < n_nodes; i ++)
	{
		exit_if_false(node[i] = node_new(i),"allocating a node");
		exit_if_false(node_read_x(input_file,node[i]) == NODE_SUCCESS,"reading a node position");
	}

	// read the faces
	exit_if_false(fsetpos(input_file,&position[1]) == 0,"finding face data");
	exit_if_false(fscanf(input_file,"%i",&n_faces) == 1,"reading the number of faces");
	exit_if_false(face = (FACE *)malloc(n_faces * sizeof(FACE)),"allocating faces");
	for(i = 0; i < n_faces; i ++)
	{
		exit_if_false(face[i] = face_new(i),"allocating a face");
		exit_if_false(face_read_nodes(input_file,node,face[i]) == FACE_SUCCESS,"reading a face nodes");
	}

	// read the elements
	exit_if_false(fsetpos(input_file,&position[2]) == 0,"finding element data");
	exit_if_false(fscanf(input_file,"%i",&n_elements) == 1,"reading the number of elements");
	exit_if_false(element = (ELEMENT *)malloc(n_elements * sizeof(ELEMENT)),"allocating elements");
	for(i = 0; i < n_elements; i ++)
	{
		exit_if_false(element[i] = element_new(i),"allocating an element");
		exit_if_false(element_read_faces(input_file,face,element[i]) == ELEMENT_SUCCESS,"reading an element faces");
	}

	// read the boundaries
	exit_if_false(fsetpos(input_file,&position[3]) == 0,"finding boundary data");
	exit_if_false(fscanf(input_file,"%i",&n_boundaries) == 1,"reading the number of boundaries");
	exit_if_false(boundary = (BOUNDARY *)malloc(n_boundaries * sizeof(BOUNDARY)),"allocating boundaries");
	for(i = 0; i < n_boundaries; i ++)
	{
		exit_if_false(boundary[i] = boundary_new(i),"allocating a boundary");
		exit_if_false(boundary_read_name(input_file,boundary[i]) == BOUNDARY_SUCCESS,"reading a boundary name");
		exit_if_false(boundary_read_faces(input_file,face,boundary[i]) == BOUNDARY_SUCCESS,"reading a boundary faces");
		exit_if_false(boundary_read_condition(input_file,boundary[i]) == BOUNDARY_SUCCESS,"reading a boundary condition");
	}

	// read the constants
	exit_if_false(fsetpos(input_file,&position[4]) == 0,"finding constant data");
	for(i = 0; i < solver_n_constants(); i ++)
	{
		exit_if_false(fscanf(input_file,"%s",string) == 1,"reading a constant name");
		exit_if_false(fscanf(input_file,"%lf",&value) == 1,"reading constant %s value",string);
		exit_if_false(solver_constant_set_value(string,value),"constant %s not recognised",string);
	}
	exit_if_false(solver_constants_set(),"constants not found");

	// read the initialisation
	// if( NO DATA FILE )
	exit_if_false(initial = (EXPRESSION *)malloc(solver_n_variables() * sizeof(EXPRESSION)),"allocating initial expressions");
	exit_if_false(fsetpos(input_file,&position[5]) == 0,"finding initial data");
	string[0] = '\0';
	exit_if_false(solver_constants_print(string),"forming constant experssion string");
	j = strlen(string);
	for(i = 0; i < solver_n_variables(); i ++)
	{
		exit_if_false(fscanf(input_file,"%s",&string[j]) == 1,"reading an input expression");
		exit_if_false(initial[i] = expression_generate(string),"generating an input expression");
	}

	// clean up
	free(string);
	fclose(input_file);

	//-------------------------------------------------------------------//
	
	// connectivity
	for(i = 0; i < n_elements; i ++) element_add_border(element[i]);

	// face geometry
	for(i = 0; i < n_faces; i ++) face_calculate_size(face[i]);
	for(i = 0; i < n_faces; i ++) exit_if_false(face_calculate_normal(face[i]) == FACE_SUCCESS,"calculating a face normal");
	for(i = 0; i < n_faces; i ++) exit_if_false(face_calculate_centre(face[i]) == FACE_SUCCESS,"calculating a face centre");

	// quadrature
	for(i = 0; i < n_faces; i ++) exit_if_false(face_calculate_quadrature(face[i]) == ELEMENT_SUCCESS,"calculating a face quadrature");
	for(i = 0; i < n_elements; i ++) exit_if_false(element_calculate_quadrature(element[i]) == ELEMENT_SUCCESS,"calculating an element quadrature");

	// element geometry
	for(i = 0; i < n_elements; i ++) element_calculate_size(element[i]);
	for(i = 0; i < n_elements; i ++) exit_if_false(element_calculate_centre(element[i]) == ELEMENT_SUCCESS,"calculating an element centre");

	// element numerics
	exit_if_false(element_interpolation_start() == ELEMENT_SUCCESS,"initialising element interpolation");
	for(i = 0; i < n_elements; i ++) exit_if_false(element_interpolation_calculate(element[i]) == ELEMENT_SUCCESS,"calculating an element interpolation");
	element_interpolation_end();

	// face numerics
	exit_if_false(face_interpolation_start() == ELEMENT_SUCCESS,"initialising face interpolation");
	for(i = 0; i < n_faces; i ++) exit_if_false(face_interpolation_calculate(face[i]) == ELEMENT_SUCCESS,"calculating a face interpolation");
	face_interpolation_end();

	// unknowns
	for(i = 0; i < n_elements; i ++) exit_if_false(element_calculate_unknowns(element[i],n_elements) == ELEMENT_SUCCESS,"calculating an element unknowns");

	// initialise system
	SPARSE system = sparse_new(NULL);
	exit_if_false(system != NULL,"allocating the system");
	for(i = 0; i < n_elements; i ++) exit_if_false(element_add_to_system(element[i],system),"adding an element to the system");
	for(i = 0; i < n_faces; i ++) exit_if_false(face_add_to_system(face[i],system),"adding a face to the system");
	exit_if_false(sparse_order_sub_matrices(system) == SPARSE_SUCCESS,"ordering the system");

	// modularise system assembly
	// // interpolation onto point value lists
	// // calculation of point residual and jacobians
	// // construction of local matrices
	// // addition into system

	//-------------------------------------------------------------------//
	
	//for(i = 0; i < n_nodes; i ++) node_print(node[i]);
	for(i = 0; i < n_faces; i ++) face_print(face[i]);
	for(i = 0; i < n_elements; i ++) element_print(element[i]);
	//for(i = 0; i < n_boundaries; i ++) boundary_print(boundary[i]);
	
	//printf("set term wxt 0\n");
	//face_plot(face[20]);
	//printf("set term wxt 1\n");
	//element_plot(element[10]);
	
	//sparse_print(system);
	//sparse_spy(system,40,40);
	
	//-------------------------------------------------------------------//

	// clean up
	solver_end();
	for(i = 0; i < n_nodes; i ++) node_free(node[i]);
	free(node);
	for(i = 0; i < n_faces; i ++) face_free(face[i]);
	free(face);
	for(i = 0; i < n_elements; i ++) element_free(element[i]);
	free(element);
	for(i = 0; i < n_boundaries; i ++) boundary_free(boundary[i]);
	free(boundary);
	for(i = 0; i < solver_n_variables(); i ++) expression_free(initial[i]);
	free(initial);
	sparse_free(system);

	return 0;
}

//////////////////////////////////////////////////////////////////
