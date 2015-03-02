#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_sf_bessel.h>
#include "pwave_math.h"
#include "pwave_io.h"

void update_field(Field *f, Droplet *drops, size_t n_drops)
{
	double S = f->S;
	//expM = e^(-1/M)
	double expM = 0.90;
	
	//variables for the loop
	double p_x,p_y;
	double delta;
	double dist;
	for(int i=0;i<f->x_max;i++)
	{
		for(int j=0;j<f->y_max;j++)
		{
			Index ind = {i,j};
			delta = 0;
			for(int m=0;m<n_drops;m++)
			{
				p_x = drops[m].p.x;
				p_y = drops[m].p.y;
				dist = sqrt(pow(p_x-i,2)+pow(p_y-j,2));
				delta = delta + S * gsl_sf_bessel_J0(dist/(f->a));
			}
			sfv(f,ind,gfv(f,ind)*expM + delta);
		}
	}
}

void update_particle(Field *f_in, Droplet *drop)
{ 
  // Too lazy to add the spring potential in a clean way, omega (defined here) gives the (non-dimensionalized) spring potential frequency.
  double Omega = 0.05;
  // Same for center point.
  Point p_center = {50,50};
	Point pg = point_grad(f_in,drop->p);
	double Q = drop->Q;
	Point np = { 2/(Q+1)*drop->p.x - Omega*Omega/(Q+1)*(drop->p.x - p_center.x)+ (Q-1)/(Q+1)*drop->p_prev.x - (f_in->a)/(Q+1)*(pg.x),
		2/(Q+1)*drop->p.y - Omega*Omega/(Q+1)*(drop->p.y - p_center.y) + (Q-1)/(Q+1)*drop->p_prev.y - (f_in->a)/(Q+1)*(pg.y) };
	drop->p_prev.x = drop->p.x;
	drop->p_prev.y = drop->p.y;
	drop->p.x = np.x;
	drop->p.y = np.y;	
}

int main() {
	int max_iterations=1000;
	double L = 100;
	double S = .3;
	double a = 2;

	printf("\nRunning main...\n");
	Field test_field;
	initialize_field(&test_field,L,L,S,a);

	char *ofname = (char *)malloc(32*sizeof(char));	
	char *pfname = (char *)malloc(32*sizeof(char));
	char *fpfname = (char *)malloc(32*sizeof(char));
	
	size_t n_drops = 1;
	Droplet drops[1];
	int dsl = 6;
	
	double s = 0.1*L-0*(dsl/2)*a;
	for(int i=0;i<dsl;i++)
		for(int j=0;j<dsl;j++)
			initialize_droplet(&drops[i*dsl+j],&(Point){s+i*a,s+j*a},3);

	sprintf(pfname,"output/particle_paths.txt");
	sprintf(fpfname,"output/fparticle_paths.txt");
	FILE *path_file = fopen(pfname,"w");
  FILE *fpath_file = fopen(fpfname,"w");

	for(int i=0;i<max_iterations;i++)
	{
		printf(">%4d/%4d.. ",i,max_iterations);
		sprintf(ofname,"output/mat_%05d",i);
		FILE *ofile = fopen(ofname,"w");
		
		// Update all the particle positions!
		for(int m = 0; m<n_drops; m++)
		{
			update_particle(&test_field,&drops[m]);
			printf("%4d,%4d ",(int)round(drops[m].p.x),(int)round(drops[m].p.y));
			fprintf(path_file," %d,%d",(int)round(drops[m].p.x),(int)round(drops[m].p.y));
			fprintf(fpath_file," %f,%f",drops[m].p.x,drops[m].p.y);
		}
		printf("\n");
		fprintf(path_file,"\n");
		fprintf(fpath_file,"\n");
	
		// Update the field based on the new particle positions.
		update_field(&test_field, &drops[0],n_drops);
		print_field(&test_field,ofile);
		fclose(ofile);
	}

	fclose(path_file);
  fclose(fpath_file);
	free(pfname);
  free(fpfname);
	free(ofname);
	destroy_field(&test_field);
}
