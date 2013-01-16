static int n_conditions = 4;

void value_empty(double *parameter, double *x, double *n, double *v)
{
}
void value_pressure(double *parameter, double *x, double *n, double *v)
{
        v[0] = parameter[0];
}
void value_velocity(double *parameter, double *x, double *n, double *v)
{
        v[0] = parameter[0];
        v[1] = parameter[1];
}
void value_normal_velocity(double *parameter, double *x, double *n, double *v)
{
        double n_hat = sqrt(n[0]*n[0]+n[1]*n[1]);
        v[0] = - parameter[0]*n[0]/n_hat;
        v[1] = - parameter[0]*n[1]/n_hat;
}

static struct s_CONDITION condition[] =
{
        { "empty" , 0 , {} , {} , 0 , &value_empty } ,
        { "pressure" , 1 , {0} , {0} , 1 , &value_pressure } ,
        { "velocity" , 2 , {1,2} , {0,0} , 2 , &value_velocity } ,
        { "normal_velocity" , 2 , {1,2} , {0,0} , 1 , &value_normal_velocity }
};

CONDITION empty = &condition[0];
