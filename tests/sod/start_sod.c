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
#include "timestep.h"
#include "refinement.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "iterators.h"
#include "auxiliary.h"

void run_output() {
	char filename[256];
	FILE *output;
	int coords[nDim];
	float pos[nDim];
	int icell;

	/* output 1-d slice through middle of volume */
	sprintf(filename, "%s/%s_%04u.dat", output_directory, jobname, step );
	
	output = fopen( filename, "w" );
	if ( output == NULL ) {
		cart_error("Unable to open %s for writing!", filename );
	}

	coords[1] = coords[2] = num_grid/2;

	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		icell = root_cell_location( sfc_index(coords) );

		if ( icell != NULL_OCT ) {
			cell_position( icell, pos );

			fprintf( output, "%e %e %e %e\n", pos[0], cell_gas_density(icell),
				cell_momentum(icell,0), cell_gas_pressure(icell) );
		}
	}

	fclose(output);
}

void init_run() {
	int i;
	float pos[nDim];
	int level, icell;
	int num_level_cells;
	int *level_cells;

	build_cell_buffer();
	repair_neighbors();

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_position(icell, pos);
			
			cell_momentum(icell,0) = 0.0;
			cell_momentum(icell,1) = 0.0;
			cell_momentum(icell,2) = 0.0;

			cell_gas_gamma(icell) = gamma;

			if ( pos[0] < (float)num_grid/2.0 ) {
				cell_gas_density(icell) = 1.0;
				cell_gas_pressure(icell) = 1.0;
			} else {
				cell_gas_density(icell) = 0.125;
				cell_gas_pressure(icell) = 0.1;
			}

			cell_gas_internal_energy(icell) = cell_gas_pressure(icell)/( (cell_gas_gamma(icell)-1.0) );
			cell_gas_energy(icell) = cell_gas_internal_energy(icell);
		}

		cart_free( level_cells );

		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
		modify( level, 1 );
	}

	/* set time variables */
	tl[min_level] = t_init;
	dtl[min_level] = 0.0;
        choose_timestep( &dtl[min_level] );
        aexp[min_level] = 1.0;

	for ( level = min_level+1; level <= max_level; level++ ) {
		dtl[level] = 0.5*dtl[level-1];
		tl[level] = tl[min_level];
		aexp[level] = aexp[min_level];
	}
}
