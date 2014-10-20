#include <stdio.h>
#include <stdlib.h>
#include "pwave_math.h"




void print_field(Field *f_in, FILE *fp) 
{
	for(int i=0;i<f_in->x_max;i++)
	{
		for(int j=0;j<f_in->y_max;j++)
		{
			Index ind = {i,j};
			fprintf(fp,"%f ",gfv(f_in,ind));
		}
		fprintf(fp,"\n");
	}
}
