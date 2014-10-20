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
	double length_scale;
	size_t x_max;
	size_t y_max;
	double *field_array;
} Field;

void initialize_field(Field *f_in, size_t field_x, size_t field_y);

void destroy_field(Field *f_in);

double gfv(Field *f_in, Index in);

int sfv(Field *f_in, Index in, double val);

Point index_grad(Field *f_in, Index in);

Point point_grad(Field *f_in, Point pt);

double bess0(double x);
