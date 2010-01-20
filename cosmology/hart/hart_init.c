#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "auxiliary.h"
#include "tree.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
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

#include "hart_init.h"

void hart_init() {
	int i;
	char filename[256], filename2[256];

#ifdef PARTICLES
	sprintf( filename, "%s/PMcrd.DAT", output_directory );
	sprintf( filename2, "%s/PMcrs0.DAT", output_directory );

	restart_load_balance( NULL, filename, filename2 );

	read_particles( filename, filename2, NULL, NULL, 0, NULL );
	cart_debug("read in particles");
#endif

#ifdef HYDRO
	sprintf( filename, "%s/tr_ic.dat", output_directory );
	read_hart_gas_ic(filename);
	cart_debug("read in gas");

	hydro_magic( min_level );
	hydro_eos( min_level );
#endif /* HYDRO */

	units_reset();
	units_update( min_level );


	cart_debug("tl[min_level] = %f", tl[min_level] );
	cart_debug(" a[min_level] = %f", auni[min_level] );

	dtl[min_level] = 0.0;
	choose_timestep( &dtl[min_level] );

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
