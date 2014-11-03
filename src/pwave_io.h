typedef struct sim
{
	//Side length (simulation boundary)
	size_t L;
	//Number of iterations
	int N;
	//Number of droplets
	int N_d;

	//Field-specific parameter "S"
	double S;
	
	Droplet *drops;
}

void print_field(Field *f_in, FILE *fp);
