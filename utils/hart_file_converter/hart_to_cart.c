#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "defs.h"
#include "io.h"
#include "tree.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "auxiliary.h"

#include "extra/hart_io.h"

int main ( int argc, char *argv[]) {
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &local_proc_id );

	if ( argc != 3 ) {
		cart_error("Usage: mpirun -np 1 hart_to_cart input output");
	}

	strcpy( output_directory, "." );

	init_auxiliary();
	init_parallel_grid();
	init_tree();
	init_cell_buffer();

	read_hart_grid_binary( argv[1] );
	cart_debug("done reading data...");
	write_grid_binary( argv[2] );
	
	MPI_Finalize();

	return 0;
}
