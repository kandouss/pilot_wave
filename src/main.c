#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pwave_math.h"
#include "pwave_io.h"

typedef struct P
{
	double delta_t;
	double m;
	double D;
	double A;
	double K_f;
	double T_f;
	double M_e;

} P;

void update_field(Field *f_in, Point *pts, size_t n_pts)
{
	double scale=f_in->length_scale;
	// debug value for A Delta / T_f
	double C = 0.5;
	// debug value for exp(-Delta Me/T_f)
	double K = 0.9;
	// Debug value for K_f
	double Kf = 0.15;
	
	//variables for the loop
	double delta;
	double dist;
	for(int i=0;i<f_in->x_max;i++)
	{
		for(int j=0;j<f_in->y_max;j++)
		{
			Index ind = {i,j};
			delta = 0;
			for(int m=0;m<n_pts;m++)
			{
				dist = sqrt(pow(pts[m].x-scale*i,2)+pow(pts[m].y-scale*j,2));
				delta = delta + bess0(Kf * dist);
			}
			//sfv(f_in,ind,(old_value + C*(delta0+delta1+delta2))*K);
			sfv(f_in,ind,(gfv(f_in,ind) + C*delta)*K);
		}
	}
}

Point update_particle(Field *f_in, Point pt_1, Point pt_2)
{
	double M=1;
	double dt=.1;
	double g=1000;
	double D=5;
	Point pg = point_grad(f_in,pt_1);
//	printf("grad p is (%f,%f)\n",pg.x,pg.y);
	double c0 = D/dt + M/pow(dt,2);
	double c1 = 2*M/pow(dt,2);
	double c2 = D/dt - M/pow(dt,2);
	double c3 = -1*M*g;
	return (Point) { (c1*pt_1.x + c2*pt_2.x + c3*pg.x)/c0,
			 (c1*pt_1.y + c2*pt_2.y + c3*pg.y)/c0	}; 
}

int main() {
	int max_iterations=2000;
	double side_length = 1000;
	int pps = 500;

	printf("\nRunning main...\n");
	Field test_field;
	initialize_field(&test_field,pps,pps);
	test_field.length_scale = (double)side_length/pps;

	char *ofname = (char *)malloc(32*sizeof(char));
	char *infofname = (char *)malloc(32*sizeof(char));
	
	size_t n_pts = 5;
	Point pts[5] = {
		{side_length*0.52,side_length*0.52},
		{side_length*0.51,side_length*0.50},
		{side_length*0.56,side_length*0.49},
		{side_length*0.45,side_length*0.55},
		{side_length*0.47,side_length*0.47} };

	Point pts_prev[5],pts_tmp[5];

	// Zero initial momentum for the points.
	for(int m=0;m<n_pts;m++)
	{
		pts_prev[m] = (Point){pts[m].x,pts[m].y};
		pts_tmp[m] = (Point){pts[m].x,pts[m].y};
	}

	sprintf(infofname,"output/particle_paths.txt");
	FILE *info_file = fopen(infofname,"w");

	for(int i=0;i<max_iterations;i++)
	{
		printf("... %d/%d...\n",i,max_iterations);
		sprintf(ofname,"output/mat_%05d",i);
		FILE *ofile = fopen(ofname,"w");
		// Update all the particle positions!
		for(int m = 0; m<n_pts; m++)
		{
			pts_tmp[m] = pts_prev[m];
			pts_prev[m] = pts[m];
			
			pts[m] = update_particle(&test_field,pts_prev[m],pts_tmp[m]);
			fprintf(info_file,"%d,%d,",(int)(pts[m].x*pps/side_length),(int)(pts[m].y*pps/side_length));
		}
		fprintf(info_file,"\n");
	
		// Update the field based on the new particle positions.
		update_field(&test_field, &pts[0],n_pts);

		print_field(&test_field,ofile);
		fclose(ofile);
	}

	fclose(info_file);
	free(ofname);
	free(infofname);
	destroy_field(&test_field);
	return 0;
}
