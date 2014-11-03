#include "pwave_math.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

extern inline double gfv(Field *f_in, Index in);
extern inline void sfv(Field *f_in, Index in, double val);
extern inline Point index_grad(Field *f_in, Index in);
extern inline Point pp(Point p,double d);
extern inline Point ps(Point p1,Point p2);

void initialize_field(Field *f_in, size_t field_x, size_t field_y, double S_in, double a_in)
{
	f_in->x_max = field_x;
	f_in->y_max = field_y;
	f_in->S = S_in;
	f_in->a = a_in;
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

Point point_grad(Field *f_in, Point pt)
{
	double a = pt.x - (int)pt.x;
	double b = pt.y - (int)pt.y;

	Point g_nn = index_grad(f_in,(Index){(int)pt.x+1,(int)pt.y+1});
	Point g_np = index_grad(f_in,(Index){(int)pt.x+1,(int)pt.y});
	Point g_pn = index_grad(f_in,(Index){(int)pt.x,(int)pt.y+1});
	Point g_pp = index_grad(f_in,(Index){(int)pt.x,(int)pt.y});

	return ps( pp(g_nn,a*b),
		 ps( pp(g_pn,a*(1-b)),
		   ps( pp(g_np,(1-a)*b), 
			pp(g_pp,(1-a)*(1-b)) ) ) );

//	return (Point) {
//		(g_nn.x*a*b + g_pn.x*a*(1-b) + g_np.x*(1-a)*b + g_pp.x*(1-a)*(1-b)),
//		(g_nn.y*a*b + g_pn.y*a*(1-b) + g_np.y*(1-a)*b + g_pp.y*(1-a)*(1-b)) };
}
