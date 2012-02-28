#include "config.h"

#include <stdio.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro.h"
#include "io.h"
#include "io_art.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "run/hydro_step.h"

#include "hart_init.h"

void hart_init() {
	int i;
	char filename[256], filename2[256];

#ifdef PARTICLES
	sprintf( filename, "%s/PMcrd.DAT", output_directory );
	sprintf( filename2, "%s/PMcrs0.DAT", output_directory );

	restart_load_balance_art( NULL, filename, filename2 );

	read_art_particles( filename, filename2, NULL, NULL, 0, NULL );
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
