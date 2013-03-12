//////////////////////////////////////////////////////////////////

#include "isitek.h"

//////////////////////////////////////////////////////////////////

void process_geometry(int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element)
{
	int e, f, i, j;

	// allocate face borders
	for(f = 0; f < n_faces; f ++)
	{
		face[f].n_borders = MAX_FACE_N_BORDERS;
		exit_if_false(face[f].border = allocate_face_border(&face[f]),"allocating face borders");
		face[f].n_borders = 0;
	}

	// element geometry
	struct NODE * temp;
	struct FACE * g[MAX_ELEMENT_N_FACES];
	int o[MAX_ELEMENT_N_FACES];
	for(e = 0; e < n_elements; e ++)
	{
		// first face
		g[0] = element[e].face[0];
		o[0] = 0;

		// loop around generating connected list and orientations
		for(i = 1; i < element[e].n_faces; i ++)
		{
			for(j = 0; j < element[e].n_faces; j ++)
			{
				g[i] = element[e].face[j];
				if(g[i] == g[i-1]) continue;
				if(g[i-1]->node[o[i-1]] == g[i]->node[0]) { o[i] = 1; break; }
				if(g[i-1]->node[o[i-1]] == g[i]->node[1]) { o[i] = 0; break; }
			}
		}

		// signed area and centroid
		double r, a = 0.0, c[2] = {0.0,0.0};
		for(i = 0; i < element[e].n_faces; i ++)
		{
			r = g[i]->node[!o[i]]->x[0] * g[i]->node[o[i]]->x[1] - g[i]->node[o[i]]->x[0] * g[i]->node[!o[i]]->x[1];
			a += r;
			for(j = 0; j < 2; j ++) c[j] += r * ( g[i]->node[0]->x[j] + g[i]->node[1]->x[j] );
		}

		// if the area is -ve reverse all orientations
		if(a < 0) for(i = 0; i < element[e].n_faces; i ++) o[i] = !o[i];

		// convert to unsigned size and centre
		a *= 0.5;
		for(i = 0; i < 2; i ++) c[i] /= 6*a;
		for(i = 0; i < 2; i ++) element[e].centre[i] = c[i];
		element[e].size = sqrt(fabs(a));

		// add the element to the list of face borders
		for(i = 0; i < element[e].n_faces; i ++)
		{
			// if the face has no borders and this is in front then reverse the face
			if(g[i]->n_borders == 0 && o[i] == 1)
			{
				temp = g[i]->node[0];
				g[i]->node[0] = g[i]->node[1];
				g[i]->node[1] = temp;
			}

			g[i]->n_borders ++;
			g[i]->border[g[i]->n_borders-1] = &element[e];
		}

	}

	// reallocate face borders
	for(f = 0; f < n_faces; f ++)
	{
		exit_if_false(face[f].border = allocate_face_border(&face[f]),"allocating face borders");
	}

	// calculate face geometry
	double dx[2];
	for(f = 0; f < n_faces; f ++)
	{
		for(i = 0; i < 2; i ++) dx[i] = face[f].node[1]->x[i] - face[f].node[0]->x[i];

		face[f].size = sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

		// normal points from border 0 to border 1
		face[f].normal[0] = - dx[1]/face[f].size;
		face[f].normal[1] = + dx[0]/face[f].size;

		for(i = 0; i < 2; i ++) dx[i] = face[f].centre[i] = 0.5*(face[f].node[1]->x[i] + face[f].node[0]->x[i]);
	}
}

//////////////////////////////////////////////////////////////////
