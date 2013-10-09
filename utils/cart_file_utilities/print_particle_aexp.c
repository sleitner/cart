#include <stdlib.h>
#include <stdio.h>

#define nDim	3
#include "particle_io.h"

int main( int argc, char *argv[] ) {
	int i, j, k;
	int endian, nbody_flag;
	particle_header header;
	float max[2*nDim], min[2*nDim];
	FILE *input;
	int page;
	int num_pages, num_read;
	int num_parts_in_page, num_parts_per_page;
	particleid_t num_particles_total;
	float *input_page;

	if ( argc < 2 || argc > 3 ) {
		printf("Usage: print_particle_info PMcrd.DAT [PMcrs.DAT]\n");
		exit(1);
	}

	read_particle_header( argv[1], &header, &endian, &nbody_flag );

	printf("%e\n", header.aexpn );
	
	return 0;
}
