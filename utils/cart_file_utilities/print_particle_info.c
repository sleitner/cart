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
	printf("Read header, endian = %u, nbody format = %u\n", endian, nbody_flag);

	printf("aexpn    = %f\n", header.aexpn );
	printf("aexp0    = %f\n", header.aexp0 );
	printf("amplt    = %e\n", header.amplt );
	printf("astep    = %e\n", header.astep );
	printf("istep    = %u\n", header.istep );
	printf("partw    = %e\n", header.partw );
	printf("tintg    = %e\n", header.tintg );
	printf("ekin     = %e\n", header.ekin );
	printf("ekin1    = %e\n", header.ekin1 );
	printf("ekin2    = %e\n", header.ekin2 );
	printf("au0      = %e\n", header.au0 );
	printf("aeu0     = %e\n", header.aeu0 );
	printf("Nrow     = %u\n", header.Nrow );
	printf("Ngrid    = %u\n", header.Ngrid );
	printf("Nspecies = %u\n", header.Nspecies );
	printf("Nseed    = %u\n", header.Nseed );
	printf("Om0      = %f\n", header.Om0 );
	printf("Oml0     = %f\n", header.Oml0 );
	printf("hubble   = %f\n", header.hubble );
	printf("Wp5      = %f\n", header.Wp5 );
	printf("Ocurv    = %f\n", header.Ocurv );
	printf("Omb0     = %f\n", header.Omb0 );
	
	for ( i = 0; i < header.Nspecies; i++ ) {
		printf("Particle Specie %u: %u with %e mass\n", i, header.num[i], header.mass[i] );
	}

	/* read in particle data as well */
	if ( argc == 3 ) {
		num_particles_total = header.num[ header.Nspecies - 1 ];
		num_parts_per_page = header.Nrow*header.Nrow;
		num_pages = (num_particles_total-1) / num_parts_per_page + 1;

		input_page = malloc( 2*nDim*num_parts_per_page*sizeof(float) );

		/* start loading actual particle data */
		input = fopen( argv[2], "r" );
		if ( input == NULL ) {
			fprintf( stderr, "Unable to open particle file %s for reading!", argv[2] );
			exit(1);
		}

		for ( i = 0; i < num_pages; i++ ) {
			if ( i == num_pages - 1 ) {
				num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
			} else {
				num_parts_in_page = num_parts_per_page;
			}

			num_read = fread( input_page, sizeof(float), 2*nDim*num_parts_per_page, input );
			if ( num_read != 2*nDim*num_parts_per_page ) {
				fprintf(stderr,
					"Error reading from particle file %s: insufficient data", argv[2] );
				exit(1);
			}

			if ( endian ) {
				for ( j = 0; j < num_parts_per_page; j++ ) {
					reorder( (char *)&input_page[j], sizeof(float) );
				}
			}

			for ( j = 0; j < 2*nDim; j++ ) {
				for ( k = 0; k < num_parts_in_page; k++ ) {
					if ( input_page[num_parts_per_page*j+k] > max[j] ) {
						max[j] = input_page[num_parts_per_page*j+k];
					}

					if ( input_page[num_parts_per_page*j+k] < min[j] ) {
						min[j] = input_page[num_parts_per_page*j+k];
					}
				}
			}
		}

		printf("Position min/max: (%e,%e,%e)-(%e,%e,%e)\n", min[0], min[1], min[2],
			max[0], max[1], max[2] );
		printf("Velocity min/max: (%e,%e,%e)-(%e,%e,%e)\n", min[3], min[4], min[5],
			max[3], max[4], max[5] );

		free( input_page );
		fclose(input);
	}

	return 0;
}
