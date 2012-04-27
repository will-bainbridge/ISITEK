void print_structures(int n, int f, int e, int b, int t, int n_variables, int *variable_order, int n_nodes, struct NODE *node, int n_faces, struct FACE *face, int n_elements, struct ELEMENT *element, int n_boundaries, struct BOUNDARY *boundary, int n_terms, struct TERM *term)
{
	int i, j, k, l;

	int max_variable_order = 0;
	for(i = 0; i < n_variables; i ++) max_variable_order = MAX(max_variable_order,variable_order[i]);
	int max_n_basis = ORDER_TO_N_BASIS(max_variable_order);
	int n_gauss = ORDER_TO_N_GAUSS(max_variable_order), n_hammer = ORDER_TO_N_HAMMER(max_variable_order);

	printf("\n");

	if(f >= 0)
	{
		printf("FACE %i\n\n",f);

		printf("nodes\n");
		for(i = 0; i < face[f].n_nodes; i ++) printf("%i ",(int)(face[f].node[i] - &node[0]));
		printf("\n\n");

		printf("borders\n");
		for(i = 0; i < face[f].n_borders; i ++) printf("%i ",(int)(face[f].border[i] - &element[0]));
		printf("\n\n");

		printf("boundaries\n");
		for(i = 0; i < n_variables; i ++) {
			if(!face[f].n_boundaries[i]) printf("-");
			else for(j = 0; j < face[f].n_boundaries[i]; j ++) printf("%i ",(int)(face[f].boundary[i][j] - &boundary[0]));
			printf("\n");
		} printf("\n");

		printf("normal\n%+.4e %+.4e\n\n",face[f].normal[0],face[f].normal[1]);

		printf("centre\n%+.4e %+.4e\n\n",face[f].centre[0],face[f].centre[1]);

		printf("size\n%+.4e\n\n",face[f].size);

		printf("X\n");
		for(i = 0; i < n_gauss; i ++) {
			printf("%+.4e %+.4e\n",face[f].X[0][i],face[f].X[1][i]); 
		} printf("\n\n");

		printf("W\n");
		for(i = 0; i < n_gauss; i ++) {
			printf("%+.4e\n",face[f].W[i]); 
		} printf("\n\n");

		for(i = 0; i < n_variables; i ++) {
			printf("Q%i\n",i);
			for(j = 0; j < ORDER_TO_N_BASIS(variable_order[i]); j ++) {
				for(k = 0; k < face[f].n_borders*ORDER_TO_N_BASIS(variable_order[i]) + face[f].n_boundaries[i]; k ++) {
					for(l = 0; l < n_gauss; l ++) {
						printf("%+.4e ",face[f].Q[i][j][k][l]);
					} printf("\n");
				} printf("\n");
			}
		}
	}

	if(e >= 0)
	{
		printf("ELEMENT %i\n\n",e);

		printf("faces\n");
		for(i = 0; i < element[e].n_faces; i ++) printf("%i ",(int)(element[e].face[i] - &face[0]));
		printf("\n\n");

		printf("centre\n%+.4e %+.4e\n\n",element[e].centre[0],element[e].centre[1]);

		printf("size\n%+.4e\n\n",element[e].size);

		printf("unknowns\n");
		for(i = 0; i < n_variables; i ++) {
			for(j = 0; j < ORDER_TO_N_BASIS(variable_order[i]); j ++) {
				printf("%i ",element[e].unknown[i][j]);
			} printf("\n");
		} printf("\n");

		printf("X\n");
		for(i = 0; i < (element[e].n_faces - 2)*n_hammer; i ++) {
			printf("%+.4e %+.4e\n",element[e].X[0][i],element[e].X[1][i]); 
		} printf("\n\n");

		printf("W\n");
		for(i = 0; i < (element[e].n_faces - 2)*n_hammer; i ++) {
			printf("%+.4e\n",element[e].W[i]); 
		} printf("\n\n");

		printf("P\n");
		for(i = 0; i < max_n_basis; i ++) {
			for(j = 0; j < max_n_basis; j ++) {
				for(k = 0; k < (element[e].n_faces - 2)*n_hammer; k ++) {
					printf("%+.4e ",element[e].P[i][j][k]);
				} printf("\n");
			} printf("\n");
		}

		printf("Q\n");
		for(i = 0; i < element[e].n_faces; i ++) {
			for(j = 0; j < max_n_basis; j ++) {
				for(k = 0; k < n_gauss; k ++) {
					printf("%+.4e ",element[e].Q[i][j][k]);
				} printf("\n");
			} printf("\n");
		}

		printf("I\n");
		for(i = 0; i < n_variables; i ++) {
			for(j = 0; j < (element[e].n_faces - 2)*n_hammer; j ++) {
				for(k = 0; k < ORDER_TO_N_BASIS(variable_order[i]); k ++) {
					printf("%+.4e ",element[e].I[i][j][k]);
				} printf("\n");
			} printf("\n");
		}

		printf("V\n");
		for(i = 0; i < max_n_basis; i ++) {
			for(j = 0; j < element[e].n_faces; j ++) {
				printf("%+.4e ",element[e].V[i][j]);
			} printf("\n");
		} printf("\n");

		printf("L\n");
		for(i = 0; i < n_variables; i ++) {
			for(j = 0; j < ORDER_TO_N_BASIS(variable_order[i]); j ++) {
				for(k = 0; k < ORDER_TO_N_BASIS(variable_order[i]); k ++) {
					printf("%+.4e ",element[e].L[i][j][k]);
				} printf("\n");
			} printf("\n");
		}
	}
}
