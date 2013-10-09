#include <stdlib.h>
#include <stdio.h>

#define nDim	3
#include "particle_io.h"

int main( int argc, char *argv[] ) {
	int i, j, k;
	int endian, nbody_flag;
	particle_header header;
	FILE *input, *output;
	int page;
	int num_pages, num_read;
	int num_parts_in_page, num_parts_per_page;
	particleid_t num_particles_total;
	float *input_page;

	if ( argc < 3 || argc > 5 ) {
		printf("Usage: convert_particle_endianness PMcrd.DAT PMcrd_new.DAT [PMcrs.DAT PMcrs_new.DAT]\n");
		exit(1);
	}

	read_write_particle_header( argv[1], argv[2], &header, &endian, &nbody_flag );

	/* read in particle data as well */
	if ( argc == 5 ) {
		num_particles_total = header.num[ header.Nspecies - 1 ];
		num_parts_per_page = header.Nrow*header.Nrow;
		num_pages = (num_particles_total-1) / num_parts_per_page + 1;

		input_page = malloc( 2*nDim*num_parts_per_page*sizeof(float) );

		/* start loading actual particle data */
		input = fopen( argv[3], "r" );
		if ( input == NULL ) {
			fprintf( stderr, "Unable to open particle file %s for reading!\n", argv[3] );
			exit(1);
		}

		output = fopen( argv[4], "w" );
		if ( output == NULL ) {
			fprintf( stderr, "Unable to open %s for writing!\n", argv[4] );
			exit(1);
		}

		for ( i = 0; i < num_pages; i++ ) {
			printf("page %d/%d\n", i+1, num_pages );

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
				for ( j = 0; j < 2*nDim*num_parts_per_page; j++ ) {
					reorder( (char *)&input_page[j], sizeof(float) );
				}
			}

			fwrite( input_page, sizeof(float), 2*nDim*num_parts_per_page, output );
		}

		free( input_page );
		fclose(input);
	}

	return 0;
}
