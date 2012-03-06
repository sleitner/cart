#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "sfc.h"
#include "io.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "load_balance.h"
#include "times.h"
#include "refinement.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "iterators.h"
#include "auxiliary.h"
#include "cosmology.h"

#include "../run/step.h"

extern double t_init;

void run_output() {
	int i, j;
	char filename[128];
	int icell, sfc, size;
	float value;
	float *slice;
	FILE *output;
	int coords[nDim];

	/* output a 2-d slice through the center of the box */
	sprintf( filename, "%s/%s_slice_%04u.dat", output_directory, jobname, step );
	output = fopen( filename, "w" );

	size = num_grid;
	fwrite( &size, sizeof(int), 1, output );
	size = num_grid;
	fwrite( &size, sizeof(int), 1, output );

	slice = cart_alloc(float, num_grid*num_grid );

	coords[2] = num_grid/2;
	for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
		for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
			icell = root_cell_location( sfc_index( coords ) );
			slice[ coords[1]*num_grid + coords[0] ] = cell_gas_density(icell);
		}
	}

	fwrite( slice, num_grid*num_grid, sizeof(float), output );

	cart_free( slice );
	fclose(output);
}

void init_run() {
	int i;
	double pos[nDim];
	int level, icell;
	int num_level_cells;
	int *level_cells;

	if ( num_procs > 1 ) {
		cart_error("KH is designed to be run serially!");
	}

	/* set units to cgs */
	units_set( 1.0, 1.0, 1.0 );

	/* create cell buffer */
	build_cell_buffer();
	repair_neighbors();

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_center_position(icell, pos);

			if ( fabs(pos[1]/(double)num_grid - 0.5) < 0.25 ) {
				cell_gas_density(icell) = 2.0;
				cell_momentum(icell,0) = 0.5*cell_gas_density(icell);
			} else {
				cell_gas_density(icell) = 1.0;
				cell_momentum(icell,0) = -0.5*cell_gas_density(icell);
			}

			/* add velocity perturbation from Abel (2010) */
			cell_momentum(icell,1) = 0.1*cell_gas_density(icell) * (
					sin(4.*M_PI*(pos[0]/(double)num_grid + 0.5))*exp(-(10.0*pow(pos[1]/(double)num_grid-0.25,2.0))) +
					sin(4.*M_PI*(pos[0]/(double)num_grid))*exp(-(10.0*pow((pos[1]/(double)num_grid)-0.75,2.0))));

			cell_momentum(icell,2) = 0.0;

			cell_gas_gamma(icell) = constants->gamma;
			cell_gas_pressure(icell) = 2.5;

			cell_gas_internal_energy(icell) = cell_gas_pressure(icell)/( (cell_gas_gamma(icell)-1.0) );
			cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
		}

		cart_free( level_cells );

		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
#ifdef REFINEMENT
		modify( level, 1 );
#endif
	}

#ifdef HYDRO_TRACERS
    cart_debug("setting hydro tracers");
    set_hydro_tracers( min_level );
#endif /* HYDRO_TRACERS */

	units_init();

	/* set time variables */
	tl[min_level] = t_init;
}
