#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include "defs.h"
#include "sfc.h"
#include "io.h"
#include "units.h"
#include "constants.h"
#include "auxiliary.h"
#include "cooling.h"
#include "analysis.h"

#include "halos.h"

int main ( int argc, char *argv[]) {
	int i, j, k;
	int arg;
	FILE *input;
	float aexp;
	double a,b;
	char filename[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
	char visualization_directory[256];
	char halofinder_directory[256];
	particle_header header;
	halo_list *halos;
	int cooling_enabled = 0;
	int endian, nbody_flag;
	int coords[nDim];
	int count;
	int index;
	int sx, sy, sz;
	int num_sfcs;
	int *order;
	int *sfc_list;
	int *section_count;
	int section, icell, ihalo;
	int block_size;
	int dx, dy, dz;
	int units_set;
	int hid;
	float hx, hy, hz;

	if ( argc < 6 ) {
		cart_error("Usage: ./analyze Lbox jobname output_directory halofinder_directory aexpn1 ... aexpnN");
	}

	cart_debug("my local pid = %u", getpid() );

	units_set = 0;

	/* set up mpi datatypes, timers, units, etc */
	Lbox = atof( argv[1] );
	strcpy( jobname, argv[2] );
	strcpy( output_directory, argv[3] );
	strcpy( halofinder_directory, argv[4] );
	strcpy( visualization_directory, argv[5] );
	
	init_auxiliary();

	aexp = atof( argv[7] );
	cart_debug("analyzing a = %6.4f", aexp );

#ifdef PARTICLES
	sprintf( filename,  "%s/PMcrda%s.DAT", output_directory, argv[7] );
	read_particle_header( filename, &header, &endian, &nbody_flag );

	Omega0 = header.Om0;
	Omegab0 = header.Omb0;
	OmegaL0 = header.Oml0;
	hubble = header.hubble;
	num_grid = header.Ngrid;

	aexpn = header.aexpn;

	cart_debug("Setting units from particle file %s", filename );

	init_units();
	units_set = 1;
#endif /* PARTICLES */

	halos = load_halo_sample_file(argv[6]);
	sprintf( filename, "%s/hlist_%s.dat", halofinder_directory, argv[7] );
	load_halo_finder_catalog(filename, halos);

/*
	halos = cart_alloc( sizeof(halo_list) );
	halos->num_halos = 1;
	halos->list = cart_alloc( halos->num_halos*sizeof(halo_struct) );

	// CL10 0.5268 
	halos->list[0].id = 2;
	halos->list[0].pos[0] = 8.31349/r0;
	halos->list[0].pos[1] = 27.42838/r0;
	halos->list[0].pos[2] = 55.96809/r0;
*/

	for ( i = 0; i < halos->num_halos; i++ ) {
		halos->list[i].analysis_radius = 2.0*NUM_PIXELS*PIXEL_SIZE*r0;
	}

	/* compute halo properties */
	compute_halo_properties( output_directory, visualization_directory, halos );
	cart_debug("Done computing halo properties");

	destroy_halo_list(halos);

	return 0;
}
