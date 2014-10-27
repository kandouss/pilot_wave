#include "pwave_math.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl_sf_bessel.h>

void initialize_field(Field *f_in, size_t field_x, size_t field_y, double S_in)
{
	f_in->x_max = field_x;
	f_in->y_max = field_y;
	f_in->S = S_in;
	f_in->field_array = (double *)calloc(field_x*field_y,sizeof(double));
}

void destroy_field(Field *f_in)
{
	free(f_in->field_array);
}

void initialize_droplet(Droplet *drop, Point *start_point, double Q_in)
{
	drop->p.x = start_point->x;
	drop->p.y = start_point->y;
	
	drop->p_prev.x = start_point->x;
	drop->p_prev.y = start_point->y;

	drop->Q = Q_in;
}

// gfv stands for "get field value"
double gfv(Field *f_in, Index in)
{
	if( (in.x >= f_in->x_max) || (in.y >= f_in->y_max) )
		return -1;
	else
		return f_in->field_array[in.x*f_in->y_max + in.y];
}


//sfv stands for "set field value"
int sfv(Field *f_in, Index in, double val)
{
	if( (in.x >= f_in->x_max) || (in.y >= f_in->y_max) )
		return -1;
	else
		f_in->field_array[in.x*f_in->y_max + in.y] = val;
		return 1;
}

Point index_grad(Field *f_in, Index in)
{
	return (Point){0.5*(gfv(f_in,(Index){in.x+1,in.y})-gfv(f_in,(Index){in.x-1,in.y})),
	 0.5*(gfv(f_in,(Index){in.x,in.y+1})-gfv(f_in,(Index){in.x,in.y-1}))	};
}

Point point_grad(Field *f_in, Point pt)
{
	int x_n = (int)(pt.x+1.0);
	int x_p = (int)(pt.x);

	int y_n = (int)(pt.y+1.0);
	int y_p = (int)(pt.y);

	double a = pt.x - x_p;
	double b = pt.y - y_p;

	Point g_nn = index_grad(f_in,(Index){x_n,y_n});
	Point g_np = index_grad(f_in,(Index){x_n,y_p});
	Point g_pn = index_grad(f_in,(Index){x_p,y_n});
	Point g_pp = index_grad(f_in,(Index){x_p,y_p});
//	Weights are:
//	nn	a*b	
//	np	a*(1-b)
//	pn	(1-a)*b
//	pp	(1-a)*(1-b)
	return (Point) {
		(g_nn.x*a*b + g_np.x*a*(1-b) + g_pn.x*(1-a)*b + g_pp.x*(1-a)*(1-b)),
		(g_nn.y*a*b + g_np.y*a*(1-b) + g_pn.y*(1-a)*b + g_pp.y*(1-a)*(1-b)) };
}
double bess0(double x)
{
		return gsl_sf_bessel_J0(x);
}
