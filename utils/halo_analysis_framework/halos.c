#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "auxiliary.h"
#include "units.h"
#include "constants.h"
#include "sfc.h"

#include "halos.h"

/* loads a file with a single list of integer halo ids, used to limit which halo positions are 
 * loaded for analysis */
halo_list *load_halo_sample_file( char *filename ) {
	int i, h;
	FILE *input;
	char line[256];
	int hid;
	int hpid;
	float hx, hy, hz, hvx, hvy, hvz;
	float hrvir, hrhalo, hmvir;
	float rt, rvir,  mvir, vmax, rmax, rs;
	int np, coords[nDim];
	halo_list *halos;

	halos = cart_alloc( sizeof(halo_list) );

	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s", filename );
	}

	/* count halos */
	halos->num_halos = 0;
	while (	fscanf( input, "%u\n", &hid ) != EOF ) {
		halos->num_halos++;
	}

	cart_debug("loading sample of %u halos...", halos->num_halos );

	rewind( input );

	halos->list = cart_alloc( halos->num_halos*sizeof(halo_struct) );

	h = 0;
	while ( fscanf( input, "%u\n", &hid ) != EOF ) {
		halos->list[h].id = hid;
		h++;
	}

	cart_assert( h == halos->num_halos );

	fclose(input);

	return halos;
}

void load_halo_finder_catalog( char *filename, halo_list *halos ) {
	int i, h;
	FILE *input;
	char line[256];
	int hid;
	int hpid;
	float hx, hy, hz, hvx, hvy, hvz;
	float hrvir, hrhalo, hmvir;
	float rt, rvir,  mvir, vmax, rmax, rs;
	int np, coords[nDim];
	int found;

	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s", filename );
	}

	if ( halos->num_halos == 0 ) {
		/* want to load entire set of halos */
		/* skip header lines */
		for ( i = 0; i < 17; i++ ) {
			fgets( line, 1024, input );
		}

		while ( fscanf( input, "%u %e %e %e %e %e %e %e %e %e %u %e %e %e %u",
					&hid, &hx, &hy, &hz, &hvx, &hvy, &hvz, &hrvir, &hrhalo,
					&hmvir, &np, &vmax, &rmax, &rs, &hpid ) != EOF ) {
			halos->num_halos++;
		}

		rewind(input);

		halos->list = cart_alloc( halos->num_halos*sizeof(halo_struct) );

		for ( i = 0; i < 17; i++ ) {
			fgets( line, 1024, input );
		}

		h = 0;
		while ( fscanf( input, "%u %e %e %e %e %e %e %e %e %e %u %e %e %e %u",
					&hid, &hx, &hy, &hz, &hvx, &hvy, &hvz, &hrvir, &hrhalo,
					&hmvir, &np, &vmax, &rmax, &rs, &hpid ) != EOF ) {
			halos->list[h++].id = hid;
		}
	}
	
	cart_debug("loading %u halo positions and velocities...", halos->num_halos );

	/* skip h ader lines */
	for ( i = 0; i < 17; i++ ) {
		fgets( line, 1024, input );
	}

	h = 0;
	while ( fscanf( input, "%u %e %e %e %e %e %e %e %e %e %u %e %e %e %u",
			&hid, &hx, &hy, &hz, &hvx, &hvy, &hvz, &hrvir, &hrhalo,
			&hmvir, &np, &vmax, &rmax, &rs, &hpid ) != EOF ) {

		/* n^2 in the number of halos, but doesn't assume halo list is ordered */
		for ( i = 0; i < halos->num_halos; i++ ) {
			if ( halos->list[i].id == hid ) {
				/* convert position to code units */
				halos->list[i].pos[0] = hx/r0;
				halos->list[i].pos[1] = hy/r0;
				halos->list[i].pos[2] = hz/r0;

				halos->list[i].vel[0] = hvx/v0;
				halos->list[i].vel[1] = hvy/v0;
				halos->list[i].vel[2] = hvz/v0;

				halos->list[i].rvir = hrvir/( 1e3*r0 );
				halos->list[i].rhalo = hrhalo/( 1e3*r0 );
				halos->list[i].mvir = hmvir;

				cart_debug("halo[%u] @ %e %e %e", halos->list[i].id, halos->list[i].pos[0],
					halos->list[i].pos[1], halos->list[i].pos[2] );

				break;
			}
		}
	}

	fclose(input);
}

/* loads rvir and mvir computed using full baryonic mass */
void load_baryon_catalog( char *filename, halo_list *halos ) {
	int i, h;
	FILE *input;
	char line[1024];
	int hid;
	float hrvir, hrhalo, hmvir;
	float rt, rvir,  mvir, vmax, rmax, rs;

	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s", filename );
	} 

	while ( fgets( line, 1024, input ) ) {
		if ( line[0] != '#' ) {
			sscanf( line, "%u %e %*e %*e %*e %*e %*e %*e %*e %e",
					&hid, &hrvir, &hmvir );

			for ( h = 0; h < halos->num_halos; h++ ) {
				if ( halos->list[h].id == hid ) {
					halos->list[h].rvir = hrvir/r0/1000.0;
					halos->list[h].mvir = hmvir;
					break;
				}
			}
		}
	}

	fclose(input);
}

void destroy_halo_list( halo_list *halos ) {
	if ( halos->list != NULL ) {
		cart_free( halos->list );
	}
	cart_free( halos );
}
