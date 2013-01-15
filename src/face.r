struct s_FACE
{
	int index;

	NODE *node;

	int n_borders;
	ELEMENT *border;

	BOUNDARY boundary;

	double *normal;
	double *centre;
	double size;

	int n_quadrature;
	double **X; // 2 * n_quadrature
	double *W; // n_quadrature

	int system_index;

	double ***Q; // n_interpolations * ( n_borders * n_bases + condition_n_variables ) * n_quadrature
};
