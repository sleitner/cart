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
#include "analysis_xray.h"

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
	halo_list *halos, *subhalos;
	int xray_enabled = 0;
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

	if ( argc != 8) {
		cart_error("Usage: ./analyze Lbox jobname output_directory halofinder_directory analysis_directory samplefile aexpn");
	}

	cooling_enabled = 0;
	xray_enabled = 0;
	units_set = 0;

	/* set up mpi datatypes, timers, units, etc */
	Lbox = atof( argv[1] );
	strcpy( jobname, argv[2] );
	strcpy( output_directory, argv[3] );
	strcpy( halofinder_directory, argv[4] );
	strcpy( analysis_directory, argv[5] );

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
	init_sfc();

	units_set = 1;
#else 
	init_sfc();
#endif /* PARTICLES */

#ifdef COOLING
	if ( !cooling_enabled ) {
		init_cooling();
		cooling_enabled = 1;
	}

	set_cooling_redshift( aexpn );
#endif /* COOLING */

#ifdef ANALYSIS_XRAY
	if ( !xray_enabled ) {
		init_xray_tables();
		xray_enabled = 1;
	}

	cart_debug("set_xray_redshift aexpn = %e", aexpn );
	set_xray_redshift( aexpn );
#endif /* ANALYZE_XRAY */

	halos = load_halo_sample_file(argv[6]);
	sprintf( filename, "%s/hlist_%s.dat", halofinder_directory, argv[7] );
	load_halo_finder_catalog(filename, halos);	

	subhalos = cart_alloc( sizeof(halo_list) );
	subhalos->num_halos = 0;

#ifdef EXCLUDE_SUBHALOS
	load_baryon_catalog( filename, subhalos );
#endif
	
	/* compute halo properties */
	compute_halo_properties( output_directory, analysis_directory, halos, subhalos );
	cart_debug("done computing halo properties");

	destroy_halo_list(halos);

	return 0;
}
