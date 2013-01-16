static int n_variables = 3;
static char variable_name[][128] = { "pressure" , "x-velocity" , "y-velocity" };
static int variable_order[] = { 2 , 3 , 3 };

static int n_interpolations = 7;
static int interpolation_variable[] = { 0 , 1 , 2 , 1 , 1 , 2 , 2 };
static int interpolation_differential[] = { 0 , 0 , 0 , 1 , 2 , 1 , 2 };
static char interpolation_method[] = "iiiiiii";

static int n_constants = 2;
static char constant_name[][128] = { "density" , "viscosity" };
static double constant_value[] = { 0 , 0 };
static int constant_set[] = { 0 , 0 };

#define P u[0]
#define U u[1]
#define V u[2]
#define DU_DX u[3]
#define DU_DY u[4]
#define DV_DX u[5]
#define DV_DY u[6]

#define RHO constant_value[0]
#define MU constant_value[1]

void residual(double *u, double **r) // r[equation][direction]
{
	r[0][0] = RHO*U;
	r[0][1] = RHO*V;

	r[1][0] = P + RHO*U*U - 2*MU*DU_DX;
	r[1][1] = RHO*U*V - MU*(DU_DY+DV_DX);

	r[2][0] = RHO*U*V - MU*(DU_DY+DV_DX);
	r[2][1] = P + RHO*V*V - 2*MU*DV_DY;
}

void jacobian(double *u, double ***j) // j[equation][direction][interpolation]
{
	int i;
	for(i = 0; i < n_variables*2*n_interpolations; i ++) j[0][0][i] = 0.0;
	
	j[0][0][1] = RHO;
	j[0][1][2] = RHO;

	j[1][0][0] = 1.0;
	j[1][0][1] = 2*RHO*U;
	j[1][0][3] = - 2*MU;
	j[1][1][1] = RHO*V;
	j[1][1][2] = RHO*U;
	j[1][1][4] = - MU;
	j[1][1][5] = - MU;

	j[2][0][1] = RHO*V;
	j[2][0][2] = RHO*U;
	j[2][0][4] = - MU;
	j[2][0][5] = - MU;
	j[2][1][0] = 1.0;
	j[2][1][2] = 2*RHO*V;
	j[2][1][6] = - 2*MU;
}
