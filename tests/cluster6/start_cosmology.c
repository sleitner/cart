#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include "auxiliary.h"
#include "tree.h"
#include "particle.h"                                                                                                                              
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "times.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "gravity.h"
#include "density.h"
#include "starformation.h"
#include "io.h"

#ifdef HYDRO
void read_hart_gas_ic( char *filename );
#endif /* HYDRO */

void run_output() {
}

void init_run() {
	int i,j;
	int level;
	int total_cells_per_level[max_level-min_level+1];
	char filename[256], filename2[256];

#ifdef PARTICLES
	sprintf( filename, "ICs/PMcrd.DAT" );
	sprintf( filename2, "ICs/PMcrs0.DAT" );

	restart_load_balance( NULL, filename, filename2 );

	read_particles( filename, filename2, NULL, NULL, 0, NULL );
	cart_debug("read in particles");
#endif

#ifdef HYDRO
	sprintf( filename, "ICs/tr_ic.dat" );
	read_hart_gas_ic(filename);
	cart_debug("read in gas");

	hydro_magic( min_level );
	hydro_eos( min_level );
#endif /* HYDRO */

    units_reset();
    units_update( min_level );

	cart_debug("tl[min_level] = %f", tl[min_level] );
	cart_debug(" a[min_level] = %f", auni[min_level] );

#ifdef PARTICLES
	build_mesh();

#ifdef STARFORM
	for ( i = 0; i < nDim; i++ ) {
		star_formation_volume_min[i] = refinement_volume_min[i];
		star_formation_volume_max[i] = refinement_volume_max[i];
	}
#endif /* STARFORM */

#else
	for ( i = 0; i < nDim; i++ ) {
		refinement_volume_min[i] = 0;
		refinement_volume_max[i] = (double)num_grid;
	}
#endif /* PARTICLES */

	if ( !buffer_enabled ) {
		cart_debug("building cell buffer");
		build_cell_buffer();
		repair_neighbors();
	}     
}
