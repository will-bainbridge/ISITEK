//////////////////////////////////////////////////////////////////

#include "isitek.h"

//////////////////////////////////////////////////////////////////

void process_geometry(int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element)
{
	int e, f, i, j;

	// generate face borders
	for(f = 0; f < n_faces; f ++)
	{
		face[f].n_borders = MAX_FACE_N_BORDERS;
		exit_if_false(face[f].border = allocate_face_border(&face[f]),"allocating face borders");
		face[f].n_borders = 0;
	}
	for(e = 0; e < n_elements; e ++)
	{
		for(i = 0; i < element[e].n_faces; i ++)
		{
			element[e].face[i]->border[element[e].face[i]->n_borders ++] = &element[e];
		}
	}
	for(f = 0; f < n_faces; f ++)
	{
		exit_if_false(face[f].border = allocate_face_border(&face[f]),"allocating face borders");
	}

	// ascertain element face orientations
	int index;
	struct NODE *vertex;
	for(e = 0; e < n_elements; e ++)
	{
		exit_if_false(element[e].orient = allocate_element_orient(&element[e]),"allocating the element face orientations");

		index = 0;
		element[e].orient[0] = 0;

		for(i = 1; i < element[e].n_faces; i ++)
		{
			vertex = element[e].face[index]->node[!element[e].orient[index]];

			for(j = 0; j < element[e].n_faces; j ++)
			{
				if(j == index) continue;
				if(element[e].face[j]->node[0] == vertex || element[e].face[j]->node[1] == vertex)
				{
					index = j;
					element[e].orient[j] = element[e].face[j]->node[1] == vertex;
					break;
				}
			}

			exit_if_false(j < element[e].n_faces,"finding the next face");
		}

		for(i = 0; i < element[e].n_faces; i ++) element[e].orient[i] = 2*element[e].orient[i] - 1;
	}

	// calculate element sizes and centroids
	double area, cross;
	for(e = 0; e < n_elements; e ++)
	{
		area = 0.0;
		for(i = 0; i < 2; i ++) element[e].centre[i] = 0.0;

		for(i = 0; i < element[e].n_faces; i ++)
		{
			cross = element[e].orient[i] * (
					element[e].face[i]->node[0]->x[0] * element[e].face[i]->node[1]->x[1] -
					element[e].face[i]->node[1]->x[0] * element[e].face[i]->node[0]->x[1] );
			area += cross;
			for(j = 0; j < 2; j ++) element[e].centre[j] += cross * ( element[e].face[i]->node[0]->x[j] + element[e].face[i]->node[1]->x[j] );
		}

		area *= 0.5;
		for(i = 0; i < 2; i ++) element[e].centre[i] /= 6*area;

		if(area < 0)
		{
			area = - area;
			for(i = 0; i < element[e].n_faces; i ++) element[e].orient[i] = - element[e].orient[i];
		}

		element[e].size = sqrt(area);
	}

	// calculate face normals, centres and sizes
	double dx[2];
	for(f = 0; f < n_faces; f ++)
	{
		for(i = 0; i < 2; i ++) dx[i] = face[f].node[1]->x[i] - face[f].node[0]->x[i];

		face[f].size = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		face[f].normal[0] = + dx[1]/face[f].size;
		face[f].normal[1] = - dx[0]/face[f].size;

		for(i = 0; i < 2; i ++) dx[i] = face[f].centre[i] = 0.5*(face[f].node[1]->x[i] + face[f].node[0]->x[i]);

		face[f].size *= 0.5;
	}
}

//////////////////////////////////////////////////////////////////
