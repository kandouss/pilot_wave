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

inline double gfv(Field *f_in, Index in) { return f_in->field_array[in.x*f_in->y_max + in.y]; }

inline void sfv(Field *f_in, Index in, double val) {f_in->field_array[in.x*f_in->y_max + in.y] = val;}

inline Point ps(Point p1, Point p2){return (Point){p1.x+p2.x,p1.y+p2.y};}

inline Point pp(Point p, double d){return (Point){p.x*d,p.y*d};}

inline Point index_grad(Field *f_in, Index in) {  return (Point){
	0.5*(gfv(f_in,(Index){in.x+1,in.y})-gfv(f_in,(Index){in.x-1,in.y})),
        0.5*(gfv(f_in,(Index){in.x,in.y+1})-gfv(f_in,(Index){in.x,in.y-1})) }; }

Point point_grad(Field *f_in, Point pt);
