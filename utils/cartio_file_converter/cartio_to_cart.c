#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "defs.h"
#include "config.h"
#include "io.h"
#include "io_art.h"
#include "io_cartio.h"
#include "tree.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "timing.h"
#include "units.h"
#include "auxiliary.h"
#include "cosmology.h"
#include "particle.h"
#include "times.h"

void config_init();
void write_cartio_grid( char *filename, int num_files );

int main_analysis ( int argc, char *argv[]) {
	int index;
	int level;

	if ( num_procs > 1 ) {
		cart_error("cartio_to_cart designed for serial execution!\n");
	}

	if ( argc != 3 ) {
		cart_error("Usage: mpirun -np 1 cartio_to_cart input output[:num_output_files]");
	}

	read_cartio_restart( argv[1] );

	/* not sure this is necessary */
	units_reset();
	units_update(min_level);

#ifdef COSMOLOGY
	abox[min_level] = abox_from_tcode(tl[min_level]);
	auni[min_level] = auni_from_tcode(tl[min_level]);
#else
    abox[min_level] = auni[min_level];
#endif

	for ( level = min_level; level <= max_level; level++ ) {
		cart_debug("num_cells_per_level[%u] = %d", level, num_cells_per_level[level] );
	}

	for ( index = 0; index < strlen( argv[2] ); index++ ) {
		if ( argv[2][index] == ':' ) {
			argv[2][index] = '\0';
			num_art_output_files = atoi( &argv[2][index+1] );
			break;
        } else {
			num_art_output_files = 1;
		}
    }

	write_art_restart( );

	MPI_Finalize();

	return 0;
}
