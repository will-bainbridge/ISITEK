#define X_GT_EPS(x) (fabs(x) > 1e-8 ? (x) : 0)

void face_print(int f)
{
	int i, j, k;

	int max_variable_order = 0;
	for(i = 0; i < n_variables; i ++) max_variable_order = MAX(max_variable_order,variable_order[i]);
	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order);

	printf("face %i\n",f);

	printf("    face->node\n       ");
	for(i = 0; i < face[f].n_nodes; i ++) printf(" %i",(int)(face[f].node[i] - &node[0]));
	printf("\n");

	printf("    face->border\n       ");
	for(i = 0; i < face[f].n_borders; i ++) printf(" %i",(int)(face[f].border[i] - &element[0]));
	printf("\n");

	printf("    face->normal\n       ");
	for(i = 0; i < 2; i ++) printf(" %g",X_GT_EPS(face[f].normal[i]));
	printf("\n");

	printf("    face->centre\n       ");
	for(i = 0; i < 2; i ++) printf(" %g",face[f].centre[i]);
	printf("\n");

	printf("    face->size\n        %g\n",face[f].size);

	printf("    face->X");
	for(i = 0; i < n_gauss; i ++) { printf("\n       "); for(j = 0; j < 2; j ++) { printf(" %+e",face[f].X[j][i]); } }
	printf("\n");
	printf("    face->W");
	for(i = 0; i < n_gauss; i ++) printf("\n        %+e",face[f].W[i]);
	printf("\n");

	printf("    face->system_index\n        %i\n",n_elements+f);

	int iv[] = { 0 , 1 , 2 , 1 , 1 , 2 , 2 };
	int id[] = { 0 , 0 , 0 , 1 , 2 , 1 , 2 };

	printf("    face->Q");
	for(i = 0; i < 7; i ++) {
		for(j = 0; j < face[f].n_borders*ORDER_TO_N_BASIS(variable_order[iv[i]]) + face[f].n_boundaries[iv[i]]*n_gauss; j ++) {
			printf("\n       ");
			for(k = 0; k < n_gauss; k ++) {
				printf(" %+e",X_GT_EPS(face[f].Q[iv[i]][id[i]][j][k]));
			}
		}
		printf("\n");
	}
}

void element_print(int e)
{
	int i, j, k;

	int max_variable_order = 0;
	for(i = 0; i < n_variables; i ++) max_variable_order = MAX(max_variable_order,variable_order[i]);
	int max_n_basis = ORDER_TO_N_BASIS(max_variable_order);
	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order);

	printf("element %i\n",e);

	printf("    element->face\n       ");
	for(i = 0; i < element[e].n_faces; i ++) printf(" %i",(int)(element[e].face[i] - &face[0]));
	printf("\n");

	printf("    element->X");
	for(i = 0; i < (element[e].n_faces - 2)*n_hammer; i ++) { printf("\n       "); for(j = 0; j < 2; j ++) printf(" %+e",element[e].X[j][i]); }
	printf("\n");
	printf("    element->W");
	for(i = 0; i < (element[e].n_faces - 2)*n_hammer; i ++) printf("\n        %+e",element[e].W[i]);
	printf("\n");

	printf("    element->size\n        %g\n",element[e].size);

	printf("    element->centre\n       ");
	for(i = 0; i < 2; i ++) printf(" %g",element[e].centre[i]);
	printf("\n");

	printf("    element->unknown");
	for(i = 0; i < n_variables; i ++) { printf("\n       "); for(j = 0; j < ORDER_TO_N_BASIS(variable_order[i]); j ++) printf(" %i",element[e].unknown[i][j]); }
	printf("\n");

	printf("    element->system_index\n        %i\n",e);

	printf("    element->P");
	for(i = 0; i < max_n_basis; i ++) {
		for(j = 0; j < max_n_basis; j ++) {
			printf("\n       ");
			for(k = 0; k < (element[e].n_faces - 2)*n_hammer; k ++) {
				printf(" %+e",X_GT_EPS(element[e].P[i][j][k]));
			}
		}
		printf("\n");
	}

	printf("    element->Q");
	for(i = 0; i < element[e].n_faces; i ++) {
		for(j = 0; j < max_n_basis; j ++) {
			printf("\n       ");
			for(k = 0; k < n_gauss; k ++) {
				printf(" %+e",X_GT_EPS(element[e].Q[i][j][k]));
			}
		}
		printf("\n");
	}

	printf("    element->I");
	for(i = 0; i < n_variables; i ++) {
		for(j = 0; j < (element[e].n_faces - 2)*n_hammer; j ++) {
			printf("\n       ");
			for(k = 0; k < ORDER_TO_N_BASIS(variable_order[i]); k ++) {
				printf(" %+e",X_GT_EPS(element[e].I[i][j][k]));
			}
		}
		printf("\n");
	}

	printf("    element->V");
	for(i = 0; i < max_n_basis; i ++) {
		printf("\n       ");
		for(j = 0; j < element[e].n_faces; j ++) {
			printf(" %+e",X_GT_EPS(element[e].V[i][j]));
		}
	}
	printf("\n");

	printf("    element->L");
	for(i = 0; i < n_variables; i ++) {
		for(j = 0; j < ORDER_TO_N_BASIS(variable_order[i]); j ++) {
			printf("\n       ");
			for(k = 0; k < ORDER_TO_N_BASIS(variable_order[i]); k ++) {
				printf(" %+e",X_GT_EPS(element[e].L[i][j][k]));
			}
		}
		printf("\n");
	}
}
