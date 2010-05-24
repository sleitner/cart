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
#include "analysis.h"
#include "skiplist.h"

int main ( int argc, char *argv[]) {
	int i, j, k;
	int arg;
	float aexp;
	double a,b;
	char filename[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
	char analysis_directory[256];
	char halofinder_directory[256];
	particle_header header;
	halo_list *halos;
	int units_set;
	int endian, nbody_flag;
	int coords[nDim];
	int count;
	int index;
	double radius_factor;
	char radius_name[256];

	if ( argc < 10 ) {
		cart_error("Usage: ./analyze Lbox jobname output_directory halofinder_directory analysis_directory projection_extent radius sample_list aexpn");
	}

	units_set = 0;

	/* set up mpi datatypes, timers, units, etc */
	Lbox = atof( argv[1] );
	strcpy( jobname, argv[2] );
	strcpy( output_directory, argv[3] );
	strcpy( halofinder_directory, argv[4] );
	strcpy( analysis_directory, argv[5] );

	init_auxiliary();

	aexp = atof( argv[9] );
	cart_debug("analyzing a = %6.4f", aexp );

#ifdef PARTICLES
	sprintf( filename,  "%s/PMcrda%s.DAT", output_directory, argv[9] );
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

	init_sfc();

	projection_extent = atof(argv[6]);

	/* load list of halos (can't be done until units set from read_restart) */
	cart_debug("loading halo catalog");
	sprintf( filename, "%s/hlist_%s.dat", halofinder_directory, argv[9] );
	halos = load_halo_sample_file( argv[8] );
	load_halo_finder_catalog(filename, halos );
	cart_debug("done loading halo catalog");

	cart_debug("loading halo virial radii");
	if ( !strcmp( "vir", argv[7] ) ) {
		sprintf( filename, "%s/h_blist_%s.dat", analysis_directory, argv[9] );
	} else { 
		sprintf( filename, "%s/h_blist_%s_a%s.dat", analysis_directory, 
				argv[7], argv[9] );
	}
	load_baryon_catalog( filename, halos ); 
	cart_debug("done loading halo virial radii");

	compute_halo_properties( output_directory, analysis_directory, argv[7], halos );
	cart_debug("done computing halo properties");

	destroy_halo_list(halos);
	return 0;
}
