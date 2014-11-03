#ifndef STDDEF_INCLUDED
#define STDDEF_INCLUDED
#include <stddef.h>
#endif

typedef struct Point
{
	double x;
	double y;
} Point;

typedef struct Index
{
	size_t x;
	size_t y;
} Index;

typedef struct Field
{
	//Non-dimensional field parameter
	double S;
	double a;
	//Size of the 2D field array
	size_t x_max;
	size_t y_max;
	//Pointer to field array
	double *field_array;
} Field;

typedef struct Droplet
{
	//Droplet parameter
	double Q;
	//Position of the droplet
	Point p;
	//Previous position of the droplet
	Point p_prev;
} Droplet;

void initialize_field(Field *f_in, size_t field_x, size_t field_y, double S_in, double a_in);

void destroy_field(Field *f_in);

void initialize_droplet(Droplet *drop, Point *start_point, double Q_in);

double gfv(Field *f_in, Index in);

int sfv(Field *f_in, Index in, double val);

Point index_grad(Field *f_in, Index in);

Point point_grad(Field *f_in, Point pt);

double bess0(double x);
