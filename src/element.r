struct s_ELEMENT
{
	int index;

	int n_faces;
	FACE *face;

	double *centre;
	double size;

	int n_quadrature;
	double **X; // 2 * n_quadrature
	double *W; // n_quadrature

	int **unknown; // n_variables * n_bases

	int system_index;

	double ***P; // n_differentials * n_bases * n_quadrature
	double ***Q; // n_faces * n_bases * n_quadrature

	double ***I; // n_variables * n_quadrature * n_bases

	double  **V; // n_bases * n_faces
	double ***L; // n_variables * n_bases * n_bases
};
