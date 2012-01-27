//////////////////////////////////////////////////////////////////

#include "isitek.h"

#include "constants.h"

void node_read_geometry(FILE *file, struct NODE *node);
void face_read_geometry(FILE *file, struct FACE *face, struct NODE *node);
void element_read_geometry(FILE *file, struct ELEMENT *element, struct FACE *face);

void node_write_case(FILE *file, struct NODE *node);
void node_read_case(FILE *file, struct NODE *node);
void face_write_case(FILE *file, int n_variables, int *n_basis, int n_gauss, struct NODE *node, struct FACE *face, struct ELEMENT *element, struct BOUNDARY *boundary);
void face_read_case(FILE *file, int n_variables, int *n_basis, int n_gauss, struct NODE *node, struct FACE *face, struct ELEMENT *element, struct BOUNDARY *boundary);
void element_write_case(FILE *file, int n_variables, int *n_basis, int n_gauss, int n_hammer, struct FACE *face, struct ELEMENT *element);
void element_read_case(FILE *file, int n_variables, int *n_basis, int n_gauss, int n_hammer, struct FACE *face, struct ELEMENT *element);
void boundary_write_case(FILE *file, struct FACE *face, struct BOUNDARY *boundary);
void boundary_read_case(FILE *file, struct FACE *face, struct BOUNDARY *boundary);

#define MAX_N_INDICES 100

#define MAX_N_BOUNDARIES 20
#define MAX_BOUNDARY_N_FACES 100
#define BOUNDARY_LABEL "boundary"
#define BOUNDARY_FORMAT "sisd"

#define MAX_N_TERMS 20
#define MAX_TERM_N_VARIABLES 10
#define TERM_LABEL "term"
#define TERM_FORMAT "cissssd"

//////////////////////////////////////////////////////////////////

void read_geometry(FILE *file, int *n_nodes, struct NODE **node, int *n_faces, struct FACE **face, int *n_elements, struct ELEMENT **element)
{
	char *line = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	exit_if_false(line != NULL,"allocating line string");

	int i;

	while(fgets(line, MAX_STRING_LENGTH, file) != NULL)
	{
		if(strncmp(line,"NODES",5) == 0)
		{
			exit_if_false(sscanf(&line[6],"%i",n_nodes) == 1,"reading the number of nodes");
			*node = allocate_nodes(*n_nodes);
			exit_if_false(*node != NULL,"allocating the nodes");
			for(i = 0; i < *n_nodes; i ++) node_read_geometry(file, &(*node)[i]);
		}
		if(strncmp(line,"FACES",5) == 0)
		{
			exit_if_false(sscanf(&line[6],"%i",n_faces) == 1,"reading the number of faces");
			*face = allocate_faces(*n_faces);
			exit_if_false(*face != NULL,"allocating the faces");
			for(i = 0; i < *n_faces; i ++) face_read_geometry(file, &(*face)[i], *node);
		}
		if(strncmp(line,"CELLS",5) == 0 || strncmp(line,"ELEMENTS",8) == 0)
		{
			exit_if_false(sscanf(&line[6],"%i",n_elements) == 1,"reading the number of elements");
			*element = allocate_elements(*n_elements);
			exit_if_false(*element != NULL,"allocating the elements");
			for(i = 0; i < *n_elements; i ++) element_read_geometry(file, &(*element)[i], *face);
		}
	}

	exit_if_false(*n_nodes > 0,"finding nodes in the geometry file");
	exit_if_false(*n_faces > 0,"finding faces in the geometry file");
	exit_if_false(*n_elements > 0,"finding elements in the geometry file");

	free(line);
}

//--------------------------------------------------------------//

void node_read_geometry(FILE *file, struct NODE *node)
{
	int info = fscanf(file,"%lf %lf\n",&(node->x[0]),&(node->x[1]));
	exit_if_false(info == 2, "reading a node's coordinates");
}

//--------------------------------------------------------------//

void face_read_geometry(FILE *file, struct FACE *face, struct NODE *node)
{
	int i;

	//temporary storage
	int *index, count, offset;
	char *line, *temp;
	index = (int *)malloc(MAX_FACE_N_NODES * sizeof(int));
	line = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	temp = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	exit_if_false(index != NULL && line != NULL && temp != NULL ,"allocating temporary storage");

	//read a line
	exit_if_false(fgets(line, MAX_STRING_LENGTH, file) != NULL, "reading a face line");

	//strip newlines and whitespace off the end of the line
	for(i = strlen(line)-1; i >= 0; i --) if(line[i] != ' ' && line[i] != '\n') break;
	line[i+1] = '\0';

	//sequentially read the integers on the line
	count = offset = 0;
	while(offset < strlen(line))
	{
		sscanf(&line[offset],"%s",temp);
		sscanf(temp,"%i",&index[count]);
		count ++;
		offset += strlen(temp) + 1;
		while(line[offset] == ' ') offset ++;
	}

	//number of faces
	face->n_nodes = count;

	//allocate the faces
	exit_if_false(face->node = allocate_face_node(face),"allocating face nodes");

	//node pointers
	for(i = 0; i < count; i ++) face->node[i] = &node[index[i]];

	//clean up
	free(index);
	free(line);
	free(temp);
}

//--------------------------------------------------------------//

void element_read_geometry(FILE *file, struct ELEMENT *element, struct FACE *face)
{
        int i;

        //temporary storage
        int *index, count, offset;
        char *line, *temp;
        index = (int *)malloc(MAX_ELEMENT_N_FACES * sizeof(int));
        line = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
        temp = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
        exit_if_false(index != NULL && line != NULL && temp != NULL ,"allocating temporary storage");

        //read the line
        exit_if_false(fgets(line, MAX_STRING_LENGTH, file) != NULL, "reading a element line");

        //eat up whitespace and newlines
        for(i = strlen(line)-1; i >= 0; i --) if(line[i] != ' ' && line[i] != '\n') break;
        line[i+1] = '\0';

        //sequentially read the integers on the line
        count = offset = 0;
        while(offset < strlen(line))
        {
                sscanf(&line[offset],"%s",temp);
                sscanf(temp,"%i",&index[count]);
                count ++;
                offset += strlen(temp) + 1;
                while(line[offset] == ' ') offset ++;
        }

        //number of faces
        element->n_faces = count;

        //allocate the faces
        exit_if_false(element->face = allocate_element_face(element),"allocating element faces");

        //face pointers
        for(i = 0; i < count; i ++) element->face[i] = &face[index[i]];

        //clean up
        free(index);
        free(line);
        free(temp);
}

//////////////////////////////////////////////////////////////////

void write_case(FILE *file, int n_variables, int *variable_order, int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element, int n_boundaries, struct BOUNDARY *boundary)
{
        int i;

        // number of variables and orders
        exit_if_false(fwrite(&n_variables, sizeof(int), 1, file) == 1, "writing the number of variables");
        exit_if_false(fwrite(variable_order, sizeof(int), n_variables, file) == n_variables, "writing variable orders");

	int max_variable_order = 0;
	for(i = 0; i < n_variables; i ++) max_variable_order = MAX(max_variable_order,variable_order[i]);
	int *n_basis = (int *)malloc(n_variables * sizeof(int));
	exit_if_false(n_basis != NULL,"allocating n_basis");
	for(i = 0; i < n_variables; i ++) n_basis[i] = ORDER_TO_N_BASIS(variable_order[i]);

        // numbers of structures
        exit_if_false(fwrite(&n_nodes, sizeof(int), 1, file) == 1, "writing the number of nodes");
        exit_if_false(fwrite(&n_faces, sizeof(int), 1, file) == 1, "writing the number of faces");
        exit_if_false(fwrite(&n_elements, sizeof(int), 1, file) == 1, "writing the number of elements");
        exit_if_false(fwrite(&n_boundaries, sizeof(int), 1, file) == 1, "writing the number of boundaries");

	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order);

        // structures
        for(i = 0; i < n_nodes; i ++) node_write_case(file,&node[i]);
        for(i = 0; i < n_faces; i ++) face_write_case(file,n_variables,n_basis,n_gauss,node,&face[i],element,boundary);
        for(i = 0; i < n_elements; i ++) element_write_case(file,n_variables,n_basis,n_gauss,n_hammer,face,&element[i]);
        for(i = 0; i < n_boundaries; i ++) boundary_write_case(file,face,&boundary[i]);

	free(n_basis);
}

//--------------------------------------------------------------//

void read_case(FILE *file, int *n_variables, int **variable_order, int *n_nodes, struct NODE **node, int *n_faces, struct FACE **face, int *n_elements, struct ELEMENT **element, int *n_boundaries, struct BOUNDARY **boundary)
{
        int i;

        // number of variables
        exit_if_false(fread(n_variables, sizeof(int), 1, file) == 1, "reading the number of variables");
	exit_if_false(*variable_order = (int *)realloc(*variable_order, *n_variables * sizeof(int)),"allocting variable orders");
        exit_if_false(fread(*variable_order, sizeof(int), *n_variables, file) == *n_variables, "reading variable orders");

	int max_variable_order = 0;
	for(i = 0; i < *n_variables; i ++) max_variable_order = MAX(max_variable_order,(*variable_order)[i]);
	int *n_basis = (int *)malloc(*n_variables * sizeof(int));
	exit_if_false(n_basis != NULL,"allocating n_basis");
	for(i = 0; i < *n_variables; i ++) n_basis[i] = ORDER_TO_N_BASIS((*variable_order)[i]);

        // numbers of structures
        exit_if_false(fread(n_nodes, sizeof(int), 1, file) == 1, "reading the number of nodes");
        exit_if_false(fread(n_faces, sizeof(int), 1, file) == 1, "reading the number of faces");
        exit_if_false(fread(n_elements, sizeof(int), 1, file) == 1, "reading the number of elements");
        exit_if_false(fread(n_boundaries, sizeof(int), 1, file) == 1, "reading the number of boundaries");

        exit_if_false((*node = allocate_nodes(*n_nodes)) != NULL,"allocating nodes");
        exit_if_false((*face = allocate_faces(*n_faces)) != NULL,"allocating faces");
        exit_if_false((*element = allocate_elements(*n_elements)) != NULL,"allocating elements");
        exit_if_false((*boundary = allocate_boundaries(*n_boundaries)) != NULL,"allocating boundaries");

	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order);

        // structures
        for(i = 0; i < *n_nodes; i ++) node_read_case(file,&(*node)[i]);
        for(i = 0; i < *n_faces; i ++) face_read_case(file,*n_variables,n_basis,n_gauss,*node,&(*face)[i],*element,*boundary);
        for(i = 0; i < *n_elements; i ++) element_read_case(file,*n_variables,n_basis,n_gauss,n_hammer,*face,&(*element)[i]);
        for(i = 0; i < *n_boundaries; i ++) boundary_read_case(file,*face,&(*boundary)[i]);

	free(n_basis);
}

//--------------------------------------------------------------//

void node_write_case(FILE *file, struct NODE *node)
{
	exit_if_false(fwrite(node->x, sizeof(double), 2, file) == 2, "writing the node location");
}

//--------------------------------------------------------------//

void node_read_case(FILE *file, struct NODE *node)
{
	exit_if_false(fread(node->x, sizeof(double), 2, file) == 2, "reading the node location");
}

//--------------------------------------------------------------//

void face_write_case(FILE *file, int n_variables, int *n_basis, int n_gauss, struct NODE *node, struct FACE *face, struct ELEMENT *element, struct BOUNDARY *boundary)
{
	int i, j, n;
		
	int *index = (int *)malloc(MAX_N_INDICES * sizeof(int));
	exit_if_false(index != NULL,"allocating temporary storage");

	exit_if_false(fwrite(&(face->n_nodes), sizeof(int), 1, file) == 1, "writing the number of face nodes");
	for(i = 0; i < face->n_nodes; i ++) index[i] = (int)(face->node[i] - &node[0]);
	exit_if_false(fwrite(index, sizeof(int), face->n_nodes, file) == face->n_nodes, "writing the face nodes");

	exit_if_false(fwrite(&(face->n_borders), sizeof(int), 1, file) == 1, "writing the number of face borders");
	for(i = 0; i < face->n_borders; i ++) index[i] = (int)(face->border[i] - &element[0]);
	exit_if_false(fwrite(index, sizeof(int), face->n_borders, file) == face->n_borders, "writing the face borders");

	exit_if_false(fwrite(face->n_boundaries, sizeof(int), n_variables, file) == n_variables,"writing the number of face boundaries");
	for(i = 0; i < n_variables; i ++)
	{
		if(face->n_boundaries[i])
		{
			for(j = 0; j < face->n_boundaries[i]; j ++) index[j] = (int)(face->boundary[i][j] - &boundary[0]);
			exit_if_false(fwrite(index, sizeof(int), face->n_boundaries[i], file) == face->n_boundaries[i],"writing face boundaries");
		}
	}

	exit_if_false(fwrite(face->normal, sizeof(double), 2, file) == 2, "writing the face normal");
	exit_if_false(fwrite(face->centre, sizeof(double), 2, file) == 2, "writing the face centre");
	exit_if_false(fwrite(&(face->size), sizeof(double), 1, file) == 1, "writing the face size");

	exit_if_false(fwrite(face->X[0], sizeof(double), 2*n_gauss, file) == 2*n_gauss,"writing the face integration locations");
	exit_if_false(fwrite(face->W, sizeof(double), n_gauss, file) == n_gauss,"writing the face integration weights");

	for(i = 0; i < n_variables; i ++)
	{
		n = n_basis[i]*(face->n_borders*n_basis[i] + face->n_boundaries[i])*n_gauss;
		exit_if_false(fwrite(face->Q[i][0][0], sizeof(double), n, file) == n,"writing face interpolation");
	}

	free(index);
}

//--------------------------------------------------------------//

void face_read_case(FILE *file, int n_variables, int *n_basis, int n_gauss, struct NODE *node, struct FACE *face, struct ELEMENT *element, struct BOUNDARY *boundary)
{
	int i, j, n;
		
	int *index = (int *)malloc(MAX_N_INDICES * sizeof(int));
	exit_if_false(index != NULL,"allocating temporary storage");

	exit_if_false(fread(&(face->n_nodes), sizeof(int), 1, file) == 1, "reading the number of face nodes");
	exit_if_false(face->node = allocate_face_node(face),"allocating face nodes");
	exit_if_false(fread(index, sizeof(int), face->n_nodes, file) == face->n_nodes, "reading the face nodes");
	for(i = 0; i < face->n_nodes; i ++) face->node[i] = &node[index[i]];

	exit_if_false(fread(&(face->n_borders), sizeof(int), 1, file) == 1, "reading the number of face borders");
	exit_if_false(face->border = allocate_face_border(face),"allocating face borders");
	exit_if_false(fread(index, sizeof(int), face->n_borders, file) == face->n_borders, "reading the face borders");
	for(i = 0; i < face->n_borders; i ++) face->border[i] = &element[index[i]];

	exit_if_false(face->n_boundaries = allocate_face_n_boundaries(face,n_variables),"allocating the face numbers of boundaries");
	exit_if_false(fread(face->n_boundaries, sizeof(int), n_variables, file) == n_variables,"reading the number of face boundaries");
	exit_if_false(face->boundary = allocate_face_boundary(face,n_variables),"allocating the face boundaries");
	for(i = 0; i < n_variables; i ++)
	{
		if(face->n_boundaries[i])
		{
			exit_if_false(fread(index, sizeof(int), face->n_boundaries[i], file) == face->n_boundaries[i],"reading face boundaries");
			for(j = 0; j < face->n_boundaries[i]; j ++) face->boundary[i][j] = &boundary[index[i]];
		}
	}

	exit_if_false(fread(face->normal, sizeof(double), 2, file) == 2, "reading the face normal");
	exit_if_false(fread(face->centre, sizeof(double), 2, file) == 2, "reading the face centre");
	exit_if_false(fread(&(face->size), sizeof(double), 1, file) == 1, "reading the face size");

	exit_if_false(face->X = allocate_face_x(face,n_gauss),"allocating face integration locations");
	exit_if_false(fread(face->X[0], sizeof(double), 2*n_gauss, file) == 2*n_gauss,"reading the face integration locations");
	exit_if_false(face->W = allocate_face_w(face,n_gauss),"allocating face integration weights");
	exit_if_false(fread(face->W, sizeof(double), n_gauss, file) == n_gauss,"reading the face integration weights");

	exit_if_false(face->Q = allocate_face_q(face, n_variables, n_basis, n_gauss),"allocating face interpolation");
	for(i = 0; i < n_variables; i ++)
	{
		n = n_basis[i]*(face->n_borders*n_basis[i] + face->n_boundaries[i])*n_gauss;
		exit_if_false(fread(face->Q[i][0][0], sizeof(double), n, file) == n,"reading face interpolation");
	}

	free(index);
}

//--------------------------------------------------------------//

void element_write_case(FILE *file, int n_variables, int *n_basis, int n_gauss, int n_hammer, struct FACE *face, struct ELEMENT *element)
{
	int i, n;

	int *index = (int *)malloc(MAX_N_INDICES * sizeof(int));
	exit_if_false(index != NULL,"allocating temporary storage");

	exit_if_false(fwrite(&(element->n_faces), sizeof(int), 1, file) == 1,"writing the number of element faces");
	for(i = 0; i < element->n_faces; i ++) index[i] = (int)(element->face[i] - &face[0]);
	exit_if_false(fwrite(index, sizeof(int), element->n_faces, file) == element->n_faces,"writing the element faces");
	exit_if_false(fwrite(element->orient, sizeof(int), element->n_faces, file) == element->n_faces,"writing the element orientations");

	exit_if_false(fwrite(element->centre, sizeof(double), 2, file) == 2,"writing the element centre");
	exit_if_false(fwrite(&(element->size), sizeof(double), 1, file) == 1,"writing the element size");

	for(i = 0; i < n_variables; i ++) exit_if_false(fwrite(element->unknown[i], sizeof(int), n_basis[i], file) == n_basis[i],"writing element unknowns");

	int n_points = (element->n_faces-2)*n_hammer;

	exit_if_false(fwrite(element->X[0], sizeof(double), 2*n_points, file) == 2*n_points,"writing the element integration locations");
	exit_if_false(fwrite(element->W, sizeof(double), n_points, file) == n_points,"writing the element integration weights");

	int max_n_basis = 0;
	for(i = 0; i < n_variables; i ++) max_n_basis = MAX(max_n_basis,n_basis[i]);
	
	n = max_n_basis*max_n_basis*n_points;
	exit_if_false(fwrite(element->P[0][0], sizeof(double), n, file) == n,"writing element interior interpolaton");
	n = element->n_faces*max_n_basis*n_gauss;
	exit_if_false(fwrite(element->Q[0][0], sizeof(double), n, file) == n,"writing element exterior interpolaton");

	free(index);
}

//--------------------------------------------------------------//

void element_read_case(FILE *file, int n_variables, int *n_basis, int n_gauss, int n_hammer, struct FACE *face, struct ELEMENT *element)
{
	int i, n;

	int *index = (int *)malloc(MAX_N_INDICES * sizeof(int));
	exit_if_false(index != NULL,"allocating temporary storage");

	exit_if_false(fread(&(element->n_faces), sizeof(int), 1, file) == 1,"reading the number of element faces");
	exit_if_false(element->face = allocate_element_face(element),"allocating element faces");
	exit_if_false(fread(index, sizeof(int), element->n_faces, file) == element->n_faces,"reading the element faces");
	for(i = 0; i < element->n_faces; i ++) element->face[i] = &face[index[i]];
	exit_if_false(element->orient = allocate_element_orient(element),"allocating element orientations");
	exit_if_false(fread(element->orient, sizeof(int), element->n_faces, file) == element->n_faces,"reading the element orientations");

	exit_if_false(fread(element->centre, sizeof(double), 2, file) == 2,"reading the element centre");
	exit_if_false(fread(&(element->size), sizeof(double), 1, file) == 1,"reading the element size");

	exit_if_false(element->unknown = allocate_element_unknown(element,n_variables,n_basis),"allocating element unknowns");
	for(i = 0; i < n_variables; i ++) exit_if_false(fread(element->unknown[i], sizeof(int), n_basis[i], file) == n_basis[i],"reading element unknowns");

	int n_points = (element->n_faces-2)*n_hammer;

        exit_if_false(element->X = allocate_element_x(element,n_points),"allocating element integration locations");
	exit_if_false(fread(element->X[0], sizeof(double), 2*n_points, file) == 2*n_points,"reading the element integration locations");
        exit_if_false(element->W = allocate_element_w(element,n_points),"allocating element integration weights");
	exit_if_false(fread(element->W, sizeof(double), n_points, file) == n_points,"reading the element integration weights");

	int max_n_basis = 0;
	for(i = 0; i < n_variables; i ++) max_n_basis = MAX(max_n_basis,n_basis[i]);

	n = max_n_basis*max_n_basis*n_points;
	exit_if_false(element->P = allocate_element_p(element, max_n_basis, n_points),"allocating element interior interpolaton");
	exit_if_false(fread(element->P[0][0], sizeof(double), n, file) == n,"reading element interior interpolaton");
	n = element->n_faces*max_n_basis*n_gauss;
	exit_if_false(element->Q = allocate_element_q(element, max_n_basis, n_gauss),"allocating element exterior interpolaton");
	exit_if_false(fread(element->Q[0][0], sizeof(double), n, file) == n,"reading element exterior interpolaton");

	free(index);
}

//--------------------------------------------------------------//

void boundary_write_case(FILE *file, struct FACE *face, struct BOUNDARY *boundary)
{
	int i;

	int *index = (int *)malloc(MAX_N_INDICES * sizeof(int));
	exit_if_false(index != NULL,"allocating temporary storage");

	exit_if_false(fwrite(&(boundary->n_faces), sizeof(int), 1, file) == 1,"writing the number of boundary faces");
	for(i = 0; i < boundary->n_faces; i ++) index[i] = (int)(boundary->face[i] - &face[0]);
	exit_if_false(fwrite(index, sizeof(int), boundary->n_faces, file) == boundary->n_faces,"writing the boundary faces");

	exit_if_false(fwrite(&(boundary->variable), sizeof(int), 1, file) == 1,"writing the boundary variable");
	exit_if_false(fwrite(boundary->condition, sizeof(int), 2, file) == 2,"writing the boundary condition");
	exit_if_false(fwrite(&(boundary->value), sizeof(double), 1, file) == 1,"writing the boundary value");

	free(index);
}

//--------------------------------------------------------------//

void boundary_read_case(FILE *file, struct FACE *face, struct BOUNDARY *boundary)
{
	int i;

	int *index = (int *)malloc(MAX_N_INDICES * sizeof(int));
	exit_if_false(index != NULL,"allocating temporary storage");

	exit_if_false(fread(&(boundary->n_faces), sizeof(int), 1, file) == 1,"reading the number of boundary faces");
	exit_if_false(boundary->face = allocate_boundary_face(boundary),"allocating the boundary faces");
	exit_if_false(fread(index, sizeof(int), boundary->n_faces, file) == boundary->n_faces,"reading the boundary faces");
	for(i = 0; i < boundary->n_faces; i ++) boundary->face[i] = &face[index[i]];

	exit_if_false(fread(&(boundary->variable), sizeof(int), 1, file) == 1,"reading the boundary variable");
	exit_if_false(fread(boundary->condition, sizeof(int), 2, file) == 2,"reading the boundary condition");
	exit_if_false(fread(&(boundary->value), sizeof(double), 1, file) == 1,"reading the boundary value");

	free(index);
}

//////////////////////////////////////////////////////////////////

void generate_numbered_file_path(char *file_path, char *base_path, int number)
{
	char *sub = strchr(base_path, '?');
	exit_if_false(sub != NULL,"finding substitute character \"?\" in base_path");

	*sub = '\0';
	sprintf(file_path, "%s%i%s", base_path, number, sub + 1);
	*sub = '?';
}

//////////////////////////////////////////////////////////////////

void write_data(FILE *file, int n_u, double *u, int number)
{
	exit_if_false(fwrite(&number, sizeof(int), 1, file) == 1,"writing the file number");
	exit_if_false(fwrite(&n_u, sizeof(int), 1, file) == 1,"writing the number of values");
	exit_if_false(fwrite(u, sizeof(double), n_u, file) == n_u,"writing the values");
}

//--------------------------------------------------------------//

void read_data(FILE *file, int *n_u, double **u, int *number)
{
	exit_if_false(fread(number, sizeof(int), 1, file) == 1,"reading the file number");
	exit_if_false(fread(n_u, sizeof(int), 1, file) == 1,"reading the number of values");
	exit_if_false(*u = (double *)realloc(*u, *n_u * sizeof(double)),"allocating the values");
	exit_if_false(fread(*u, sizeof(double), *n_u, file) == *n_u,"reading the values");
}

//////////////////////////////////////////////////////////////////

void boundaries_input(FILE *file, int n_faces, struct FACE *face, int *n_boundaries, struct BOUNDARY **boundary)
{
        //counters
	int i, j, n = 0, info;

	//fetch the data from the file
	FETCH fetch = fetch_new(BOUNDARY_FORMAT, MAX_N_BOUNDARIES);
	exit_if_false(fetch != NULL,"allocating boundary input");
	int n_fetch = fetch_read(file, BOUNDARY_LABEL, fetch);
	exit_if_false(n_fetch > 1,"no boundaries found in input file");
	warn_if_false(n_fetch < MAX_N_BOUNDARIES,"maximum number of boundaries reached");

	//allocate boundaries
	struct BOUNDARY *b = allocate_boundaries(n_fetch);
	exit_if_false(b != NULL,"allocating boundaries");

	//temporary storage
	int offset, index[2];
	char *range = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	char *temp = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	exit_if_false(range != NULL && temp != NULL,"allocating temporary storage");

	//consider each feteched line
	for(i = 0; i < n_fetch; i ++)
	{
		// initialise
		b[n].n_faces = MAX_BOUNDARY_N_FACES;
		b[n].face = allocate_boundary_face(&b[n]);
		exit_if_false(b->face != NULL,"allocating boundary faces");
		b[n].n_faces = 0;

		//get the range
		fetch_get(fetch, i, 0, range);

		//convert comma delimiters to whitespace
		for(j = 0; j < strlen(range); j ++) if(range[j] == ',') range[j] = ' ';

		//read the range
		offset = info = 0;
		while(offset < strlen(range))
		{
			//read the range from the string
			info = sscanf(&range[offset],"%s",temp) == 1;
			info *= sscanf(temp,"%i:%i",&index[0],&index[1]) == 2;
			warn_if_false(info,"skipping boundary with unrecognised range");
			if(!info) break;

			//store boundary in the elements in the range
			for(j = index[0]; j <= index[1]; j ++) b[n].face[b[n].n_faces ++] = &face[j];

			//move to the next range in the string
			offset += strlen(temp) + 1;
		}
		if(!info) continue;

		// re-allocate
		b[n].face = allocate_boundary_face(&b[n]);
		exit_if_false(b[n].face != NULL,"re-allocating boundary faces");

		//get the variable
		fetch_get(fetch, i, 1, &b[n].variable);

		//get the condition
		fetch_get(fetch, i, 2, temp);
		for(j = 0; j < 2; j ++) b[n].condition[j] = 0;
		if(strcmp(temp,"d") != 0)
		{
			for(j = 0; j < strlen(temp); j ++)
			{
				if(temp[j] == 'n') b[n].condition[0] ++;
				else if(temp[j] == 't') b[n].condition[1] ++;
				else { info = 0; break; }
			}
			warn_if_false(info,"skipping boundary with unrecognised condition");
		}
		if(!info) continue;

		//get the value
		fetch_get(fetch, i, 3, &b[n].value);

		//increment boundary
		n ++;
	}

	//check numbers
	warn_if_false(n == n_fetch,"skipping boundaries with unrecognised formats");

	//copy pointers
	*boundary = b;
	*n_boundaries = n;

	//clean up
	fetch_destroy(fetch);
	free(range);
	free(temp);
}

//////////////////////////////////////////////////////////////////

void terms_input(FILE *file, int *n_terms, struct TERM **term)
{
	//fetch the data
	FETCH fetch = fetch_new(TERM_FORMAT,MAX_N_TERMS);
	exit_if_false(fetch != NULL,"allocating fetch");
	int n_fetch = fetch_read(file,TERM_LABEL,fetch);
	exit_if_false(n_fetch > 1,"no terms found in input file");
	warn_if_false(n_fetch < MAX_N_TERMS,"maximum number of terms reached");

	//allocate pointers
	struct TERM *t = allocate_terms(n_fetch);
	exit_if_false(t != NULL,"allocating terms");

	//counters
	int i, j, n = 0, info;

	//temporary storage
	int var_offset, dif_offset, pow_offset, dif[2];
	char *var_string = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	char *dif_string = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	char *pow_string = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	char *temp = (char *)malloc(MAX_STRING_LENGTH * sizeof(char));
	exit_if_false(var_string != NULL && dif_string != NULL && pow_string != NULL && temp != NULL,"allocating temporary strings");
	int *vars = (int *)malloc(MAX_TERM_N_VARIABLES * sizeof(int));
	int *difs = (int *)malloc(MAX_TERM_N_VARIABLES * sizeof(int));
	int *pows = (int *)malloc(MAX_TERM_N_VARIABLES * sizeof(int));
	exit_if_false(vars != NULL && difs != NULL && pows != NULL,"allocating temporary data");

	for(i = 0; i < n_fetch; i ++)
	{
		//type
		fetch_get(fetch, i, 0, &t[n].type);
		warn_if_false(t[n].type == 's' || t[n].type == 'x' || t[n].type == 'y',"skipping term with unrecognised type");

		//equation
		fetch_get(fetch, i, 1, &t[n].equation);

		//get the variable, differential and power strings
		fetch_get(fetch, i, 2, var_string);
		fetch_get(fetch, i, 3, dif_string);
		fetch_get(fetch, i, 4, pow_string);
		for(j = 0; j < strlen(var_string); j ++) if(var_string[j] == ',') var_string[j] = ' ';
		for(j = 0; j < strlen(dif_string); j ++) if(dif_string[j] == ',') dif_string[j] = ' ';
		for(j = 0; j < strlen(pow_string); j ++) if(pow_string[j] == ',') pow_string[j] = ' ';

		//read each variable in turn
		var_offset = dif_offset = pow_offset = t[n].n_variables = 0;
		while(var_offset < strlen(var_string))
		{
			info = 1;

			//read the variable index from the string
			info *= sscanf(&var_string[var_offset],"%s",temp) == 1;
			info *= sscanf(temp,"%i",&vars[t[n].n_variables]) == 1;
			var_offset += strlen(temp) + 1;

			//read the x and y differentials and convert to a differential index
			info *= sscanf(&dif_string[dif_offset],"%s",temp) == 1;
			j = dif[0] = dif[1] = 0;
			if(info)
			{
				while(temp[j] != '\0')
				{
					dif[0] += (temp[j] == 'x');
					dif[1] += (temp[j] == 'y');
					j ++;
				}
				difs[t[n].n_variables] = powers_taylor[dif[0]][dif[1]];
			}
			dif_offset += strlen(temp) + 1;

			//read the variable powers from the string
			info *= sscanf(&pow_string[pow_offset],"%s",temp) == 1;
			info *= sscanf(temp,"%i",&pows[t[n].n_variables]) == 1;
			pow_offset += strlen(temp) + 1;

			warn_if_false(info,"skipping term with unrecognised variable format");
			if(!info) continue;

			//next variable
			t[n].n_variables ++;
		}

		//allocate the variable and differential arrays
		exit_if_false(t[n].variable = allocate_term_variable(&t[n]),"allocating term variables");
		exit_if_false(t[n].differential = allocate_term_differential(&t[n]),"allocating term differentials");
		exit_if_false(t[n].power = allocate_term_power(&t[n]),"allocating term powers");
		exit_if_false(t[n].method = allocate_term_method(&t[n]),"allocating term method");

		//copy over
		for(j = 0; j < t[n].n_variables; j ++)
		{
			t[n].variable[j] = vars[j];
			t[n].differential[j] = difs[j];
			t[n].power[j] = pows[j];
		}

		//method
		fetch_get(fetch, i, 5, t[n].method);

		//constant
		fetch_get(fetch, i, 6, &t[n].constant);

		//increment the number of terms
		n ++;
	}

	//check numbers
	warn_if_false(n_fetch == n,"skipping terms with unrecognised formats");

	//copy pointers
	*term = t;
	*n_terms = n;

	//clean up
	fetch_destroy(fetch);
	free(var_string);
	free(dif_string);
	free(pow_string);
	free(temp);
	free(vars);
	free(difs);
	free(pows);
}

//////////////////////////////////////////////////////////////////
