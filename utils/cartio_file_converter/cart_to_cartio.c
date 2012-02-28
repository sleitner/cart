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
#include "times.h"
#include "units.h"
#include "auxiliary.h"
#include "cosmology.h"
#include "particle.h"

void config_init();
void write_cartio_restart_worker( char *filename, int num_files );

int main_analysis ( int argc, char *argv[]) {
	int num_files;
	int index;
	double start, end;

#ifdef PARTICLES
#ifdef STARFORM
	 if ( argc != 7) {
        cart_error("Usage: mpirun -np X cart_to_cartio input.d[:num_input_files] PMcrd PMcrs pt stars output[:num_output_files]");
    }
#else
	if ( argc != 6 ) {
		cart_error("Usage: mpirun -np X cart_to_cartio input.d[:num_input_files] PMcrd PMcrs pt output[:num_output_files]");
	}
#endif /* STARFORM */
#else 
	if ( argc != 3 ) {
		cart_error("Usage: mpirun -np X cart_to_cartio input.d[:num_input_files] output[:num_output_files]");
	}
#endif /* PARTICLES */

	for ( index = 0; index < strlen( argv[1] ); index++ ) {
		if ( argv[1][index] == ':' ) {
			argv[1][index] = '\0';
			num_art_output_files = atoi( &argv[1][index+1] );
			break;
        } else {
			num_art_output_files = 1;
		}
    }

#ifdef PARTICLES
	restart_load_balance_art( argv[1], argv[2], argv[3] );
#else
	restart_load_balance_art(argv[1], NULL, NULL);
#endif
	read_art_grid_binary( argv[1] );

#ifdef PARTICLES
#ifdef STARFORM
	read_art_particles( argv[2], argv[3], argv[5], argv[6], 0, NULL );
#else
	read_art_particles( argv[2], argv[3], argv[4], NULL, 0, NULL );
#endif
#endif

	units_reset();
    units_update(min_level);

#ifdef COSMOLOGY
	abox[min_level] = abox_from_tcode(tl[min_level]);
	auni[min_level] = auni_from_tcode(tl[min_level]);
#else
	abox[min_level] = auni[min_level];
#endif

	cart_debug("done reading data...");

	for ( index = 0; index < strlen( argv[argc-1] ); index++ ) {
		if ( argv[argc-1][index] == ':' ) {
			argv[argc-1][index] = '\0';
			num_cartio_grid_files = atoi( &argv[argc-1][index+1] );
#ifdef PARTICLES
			num_cartio_particle_files = num_cartio_grid_files;
#endif /* PARTICLES */
			break;
		}
	}

	cart_debug("writing %d output files", num_cartio_grid_files );
	start = MPI_Wtime();
#ifdef PARTICLES
	write_cartio_restart_worker( argv[argc-1], WRITE_GRID | WRITE_PARTICLES );
#else
	write_cartio_restart_worker( argv[argc-1], WRITE_GRID );
#endif /* PARTICLES */
	end = MPI_Wtime();
	printf("%d, write cartio file, %f\n", local_proc_id, end - start);

	return 0;
}
