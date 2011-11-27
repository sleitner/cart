#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "io.h"
#include "tree.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "timing.h"
#include "units.h"
#include "auxiliary.h"
#include "cosmology.h"
#include "times.h"
#include "extra/hart_io.h"

void init_run(){}
void run_output(){}

int main ( int argc, char *argv[]) {
	MPI_Init( &argc, &argv );
	MPI_Comm_size( mpi.comm.run, &num_procs );
	MPI_Comm_rank( mpi.comm.run, &local_proc_id );

	if ( argc != 3 ) {
		cart_error("Usage: mpirun -np 1 hart_to_cart input output");
	}

	strcpy( output_directory, "." );

	init_auxiliary();
	init_timers();

	config_init();

	init_parallel_grid();
	init_tree();
	init_cell_buffer();

	read_hart_grid_binary( argv[1] );
	units_reset();
#ifdef COSMOLOGY
	abox[min_level] = abox_from_tcode(tl[min_level]);
	auni[min_level] = auni_from_tcode(tl[min_level]);
#else
	abox[min_level] = auni[min_level];
#endif
	units_update(min_level);

	cart_debug("abox = %e", abox[min_level] );
	cart_debug("done reading data...");
	cart_debug("tl = %e", tl[min_level] );

	write_grid_binary( argv[2] );
	
	MPI_Finalize();

	return 0;
}
