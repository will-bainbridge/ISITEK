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

	// face ordering and element geometry
	int index;
	double area, cross;
	struct NODE *temp;
	for(e = 0; e < n_elements; e ++)
	{
		// order the face nodes around the element
		index = 0;
		for(i = 1; i < element[e].n_faces; i ++)
		{
			for(j = 0; j < element[e].n_faces; j ++)
			{
				if(j == index) continue;
				if(element[e].face[j]->node[0] == element[e].face[index]->node[1])
				{
					index = j;
					break;
				}
				if(element[e].face[j]->node[1] == element[e].face[index]->node[1])
				{
					element[e].face[j]->node[1] = element[e].face[j]->node[0];
					element[e].face[j]->node[0] = element[e].face[index]->node[1];
					index = j;
					break;
				}
			}
		}

		// signed area and centroid
		area = 0.0;
		for(i = 0; i < 2; i ++) element[e].centre[i] = 0.0;

		for(i = 0; i < element[e].n_faces; i ++)
		{
			cross = (
					element[e].face[i]->node[0]->x[0] * element[e].face[i]->node[1]->x[1] -
					element[e].face[i]->node[1]->x[0] * element[e].face[i]->node[0]->x[1]
				);
			area += cross;
			for(j = 0; j < 2; j ++) element[e].centre[j] += cross * ( element[e].face[i]->node[0]->x[j] + element[e].face[i]->node[1]->x[j] );
		}

		area *= 0.5;
		for(i = 0; i < 2; i ++) element[e].centre[i] /= 6*area;

		// unsigned size
		element[e].size = sqrt(fabs(area));

		// reverse all faces if the area is negative
		if(area < 0)
		{
			for(i = 0; i < element[e].n_faces; i ++)
			{
				temp = element[e].face[i]->node[0];
				element[e].face[i]->node[0] = element[e].face[i]->node[1];
				element[e].face[i]->node[1] = temp;
			}
		}

		// order the face borders so that the face normal points from border 0 to border 1
		for(i = 0; i < element[e].n_faces; i ++)
		{
			if(element[e].face[i]->border[0] != &element[e])
			{
				element[e].face[i]->border[1] = element[e].face[i]->border[0];
				element[e].face[i]->border[0] = &element[e];
			}
		}
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
