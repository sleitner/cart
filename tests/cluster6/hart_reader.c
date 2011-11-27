#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "times.h"
#include "tree.h"
#include "io.h"
#include "units.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "sfc.h"
#include "iterators.h"
#include "auxiliary.h"
#include "hydro.h"
#include "starformation.h"
#include "skiplist.h"
#include "index_hash.h"
#include "refinement_indicators.h"

#ifdef HYDRO 

void read_hart_gas_ic( char *filename ) {
	int i;
	FILE *input;
	int size;
	float boxh, ainit, astep;
	int ncells;
	int endian;
	int coords[nDim];
	int index;
	int proc, icell;
	MPI_Status status;
	int page_count;
	float *input_page;
	int count[MAX_PROCS];
	float *page[MAX_PROCS];
	int *page_indices[MAX_PROCS];
	int var;
	const int num_gas_vars    = 6;
	const int var_index[] = { 	
		HVAR_GAS_DENSITY, HVAR_MOMENTUM, 
		HVAR_MOMENTUM+1, HVAR_MOMENTUM+2,
		HVAR_GAS_ENERGY, HVAR_INTERNAL_ENERGY };

	if ( local_proc_id == MASTER_NODE ) {
		input = fopen(filename, "r");
		if ( input == NULL ) {
			cart_error("Unable to open %s for reading!", filename );
		}
	
		fread( &size, sizeof(int), 1, input );
		if ( size != sizeof(float) ) {
			reorder( (char *)&size, sizeof(int) );
			endian = 1;
			if ( size != sizeof(float) ) {
				cart_error("Bad file-format in read_cell_ic");
			} 
		} else {
			endian = 0;
		}
	
		fread( &boxh, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		fread( &size, sizeof(int), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &astep, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		fread( &size, sizeof(int), 1, input );
		fread( &ncells, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&ainit, sizeof(float) );
			reorder( (char *)&astep, sizeof(float) );
			reorder( (char *)&ncells, sizeof(int) );
		}

		//Lbox = boxh;
		auni_init = ainit;

		//MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, mpi.comm.run );
		MPI_Bcast( &auni_init, 1, MPI_DOUBLE, MASTER_NODE, mpi.comm.run );
	
		cart_debug("boxh = %f", boxh );
		cart_debug("ainit = %f", ainit );
		cart_debug("astep = %f", astep );
		cart_debug("ncells = %u", ncells );

		if ( ncells != num_root_cells ) {
			cart_error("ncells in %s does not match num_root_cells (%u vs %u)", 
					filename, ncells, num_root_cells );
		}

		input_page = cart_alloc(float, num_grid*num_grid );
                for ( proc = 1; proc < num_procs; proc++ ) {
			page[proc] = cart_alloc(float, num_grid*num_grid );
			page_indices[proc] = cart_alloc(int, num_grid*num_grid );
                }

		for ( var = 0; var < num_gas_vars; var++ ) {
			for ( proc = 1; proc < num_procs; proc++ ) {
				count[proc] = 0;
			}

			fread( &size, sizeof(int), 1, input );
			for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				size = fread( input_page, sizeof(float), num_grid*num_grid, input );
				if ( size != num_grid*num_grid ) {
					cart_error("Error reading from file %s", filename );
				}

				if ( endian ) {
					for ( i = 0; i < num_grid*num_grid; i++ ) {
						reorder( (char *)&input_page[i], sizeof(float) );
					}
				}
	
				page_count = 0;
				for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
					for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
						index = sfc_index( coords );
						proc = processor_owner( index );
	
						if ( proc == local_proc_id ) {
							icell = root_cell_location( index );
							cell_var( icell, var_index[var] ) = input_page[page_count];
						} else {
							/* add cell to send page */
							page[proc][count[proc]] = input_page[page_count];
							page_indices[proc][count[proc]] = index;
							count[proc]++;
	
							if ( count[proc] == num_grid*num_grid ) {
								MPI_Send( page[proc], num_grid*num_grid, MPI_FLOAT, 
									proc, 0, mpi.comm.run );
								MPI_Send( page_indices[proc], num_grid*num_grid, 
									MPI_INT, proc, 0, mpi.comm.run );
								count[proc] = 0;
							}
						}
	
						page_count++;
					}
				}
			}
			fread( &size, sizeof(int), 1, input );
			
			/* send last variables */
			for ( proc = 1; proc < num_procs; proc++ ) {
				MPI_Send( page[proc], count[proc], MPI_FLOAT, proc, 0, mpi.comm.run );
				MPI_Send( page_indices[proc], count[proc], MPI_INT, proc, 0, mpi.comm.run );
			}
		}
	
        	fclose(input);

		cart_free( input_page );
		for ( proc = 1; proc < num_procs; proc++ ) {
			cart_free( page[proc] );
			cart_free( page_indices[proc] );
		}
	} else {
		//MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, mpi.comm.run );
		MPI_Bcast( &auni_init, 1, MPI_DOUBLE, MASTER_NODE, mpi.comm.run );

		page[local_proc_id] = cart_alloc(float, num_grid*num_grid );
		page_indices[local_proc_id] = cart_alloc(int, num_grid*num_grid );

		for ( var = 0; var < num_gas_vars; var++ ) {
			page_count = num_grid*num_grid;
			while ( page_count == num_grid*num_grid ) {
				MPI_Recv( page[local_proc_id], num_grid*num_grid, MPI_FLOAT, MASTER_NODE, 0, mpi.comm.run, &status );
				MPI_Get_count( &status, MPI_FLOAT, &page_count );
				MPI_Recv( page_indices[local_proc_id], num_grid*num_grid, MPI_INT, MASTER_NODE, 0, mpi.comm.run, &status );

				for ( i = 0; i < page_count; i++ ) {
					icell = root_cell_location( page_indices[local_proc_id][i] );
					cell_var( icell, var_index[var] ) = page[local_proc_id][i];
				}
			}
		}

		cart_free( page[local_proc_id] );
		cart_free( page_indices[local_proc_id] );
	}

	/* set gas gamma on root level */
	for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
		cell_gas_gamma(i) =  constants->gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		cell_electron_internal_energy(i) = cell_gas_internal_energy(i)*wmu/wmu_e;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
	}
}

#endif  /* HYDRO */
