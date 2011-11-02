#include "config.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "parallel.h"
#include "sfc.h"
#include "tree.h"


int num_procs;
int local_proc_id;
int proc_sfc_index[MAX_PROCS+1];


unsigned int mpi_custom_flags = MPI_CUSTOM_NONE;


void config_init_parallel()
{
  control_parameter_add2(control_parameter_int,&mpi_custom_flags,"@mpi:custom-flags","mpi_custom_flags","flags that can be set to customize MPI performance. This parameter is experimental and may be removed in the future.");
}


void config_verify_parallel()
{
}


/*******************************************************  
 * init_parallel_grid
 ******************************************************/
void init_parallel_grid() 
/* purpose: initializes proc_sfc_index array 
 * requires: local_proc_id and num_procs are already set
 */
{
	int i;
	int index;

	/* to start with, divide all root octs equally amongst all nodes */
	for ( i = 0, index = 0; i < num_procs; i++, index += num_root_cells / num_procs ) {
		proc_sfc_index[i] = index;
	}
	proc_sfc_index[num_procs] = num_root_cells;
}

/*******************************************************
 * processor_owner
 ******************************************************/
int processor_owner( int sfc ) 
/* purpose: determines which processor currently 
 * 	owns the given sfc index 
 */
{
	int a, b, c;

	cart_assert( sfc >= 0 && sfc < max_sfc_index );

	/* determine if sfc is local */
	if ( sfc < proc_sfc_index[local_proc_id] ) {
		a = 0;
		b = local_proc_id-1;
	} else if ( sfc >= proc_sfc_index[local_proc_id+1] ) {
		a = local_proc_id + 1;
		b = num_procs-1;
	} else {
		return local_proc_id;
	}

	/* do binary search between procs a & b */
	while ( a != b ) {
		c = ( a + b + 1) / 2;
		
		if ( sfc < proc_sfc_index[c] ) {
			b = c-1;
		} else {
			a = c;
		}
	}

	cart_assert( a >= 0 && a < num_procs );
	cart_assert( sfc >= proc_sfc_index[a] && sfc < proc_sfc_index[a+1] );

	return a;
}


void print_comm_contents(MPI_Comm comm, const char *name)
{
  MPI_Group world, local;
  int i, n, *ranks_local, *ranks_world;

  MPI_Comm_rank(comm,&n);

  if(n == 0)
    {
      MPI_Comm_group(MPI_COMM_WORLD,&world);
      MPI_Comm_group(comm,&local);

      MPI_Comm_size(comm,&n);

      ranks_local = cart_alloc(int,n);
      ranks_world = cart_alloc(int,n);
  
      for(i=0; i<n; i++) ranks_local[i] = i;

      MPI_Group_translate_ranks(local,n,ranks_local,world,ranks_world);

      cart_debug("Communicator %s (%p), size = %d:",name,comm,n);
      for(i=0; i<n; i++) cart_debug("id = %d -> world id = %d",i,ranks_world[i]);

      cart_free(ranks_local);
      cart_free(ranks_world);

      MPI_Group_free(&local);
      MPI_Group_free(&world);
    }
}
