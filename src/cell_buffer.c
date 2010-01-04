#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "defs.h"
#include "skiplist.h"
#include "cell_buffer.h"
#include "index_hash.h"
#include "parallel.h"
#include "sfc.h"
#include "tree.h"
#include "pack.h"
#include "iterators.h"
#include "auxiliary.h"
#include "timing.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#endif

int root_buffer_enabled = 0;
int buffer_enabled = 0;
int *buffer_cell_sfc_index;

index_hash *buffer_root_hash;
index_hash *buffer_oct_hash[MAX_PROCS];
index_hash *buffer_oct_reverse_hash[MAX_PROCS];

/* WARNING: remote_buffers[min_level] stores sfc index, while
 *  local_buffers[min_level] stores actual buffer index */
int *remote_buffers[max_level-min_level+1][MAX_PROCS];
int num_remote_buffers[max_level-min_level+1][MAX_PROCS];
int num_local_buffers[max_level-min_level+1][MAX_PROCS];
int *local_buffers[max_level-min_level+1][MAX_PROCS];

int num_buffer_cells[max_level-min_level+1];	/* number of cells we're buffering */
int buffer_oct_list[max_level-min_level+1];	/* linked list for buffered octs */

/*******************************************************
 * init_cell_buffer
 ******************************************************/
void init_cell_buffer() 
/* purpose: initializes buffer_cells array (i.e. an empty
 * 	buffer)
 */
{
	int i, j;

	buffer_enabled = 0;

	for ( i = min_level; i <= max_level; i++ ) {
		num_buffer_cells[i] = 0;
		buffer_oct_list[i] = NULL_OCT;

		for ( j = 0; j < num_procs; j++ ) {
			num_remote_buffers[i][j] = 0;
			num_local_buffers[i][j] = 0;
			remote_buffers[i][j] = NULL;
			local_buffers[i][j] = NULL;
		}
	}
}

/*******************************************************
 * cell_buffer_exists
 *******************************************************/
int cell_buffer_exists( int sfc ) 
/* purpose: determines if the specified root level cell
 *  exists in the cell buffer
 *
 *  returns: 1 if the cell exists, 0 otherwise
 */
{
	cart_assert( root_buffer_enabled );
	cart_assert( sfc >= 0 && sfc < max_sfc_index );
	return ( index_hash_lookup( buffer_root_hash, sfc ) >= 0 );
}

/*******************************************************
 * cell_buffer_local_index
 *******************************************************/
int cell_buffer_local_index( int sfc ) 
/* purpose: translates an sfc index into a local
 *  index for buffer cells
 *
 * returns: -1 if the cell is not in the buffer,
 * 	otherwise it returns the local index 
 */
{
	cart_assert( sfc >= 0 && sfc < max_sfc_index );

	if ( root_buffer_enabled ) {
		return index_hash_lookup( buffer_root_hash, sfc );
	} else {
		return -1;
	}
}

void cell_buffer_add_remote_octs( int level, int processor, int *local_octs, int num_local_octs ) {
	int i;
	int *new_buffer;

	new_buffer = cart_alloc(int, (num_remote_buffers[level][processor] + num_local_octs) );

	/* add local_octs to end of array */
	for ( i = 0; i < num_remote_buffers[level][processor]; i++ ) {
		new_buffer[i] = remote_buffers[level][processor][i];
	}

	for ( i = 0; i < num_local_octs; i++ ) {
		new_buffer[ num_remote_buffers[level][processor] + i ] = local_octs[i];
	}

	num_remote_buffers[level][processor] += num_local_octs;
	cart_free( remote_buffers[level][processor] );
	remote_buffers[level][processor] = new_buffer;
}

void cell_buffer_add_local_octs( int level, int processor, int *local_octs, int num_new_octs ) {
	int i;
	int *new_buffer;

	new_buffer = cart_alloc(int, (num_local_buffers[level][processor] + num_new_octs) );

	/* add at end (like remote octs) */
	for ( i = 0; i < num_local_buffers[level][processor]; i++ ) {
		new_buffer[i] = local_buffers[level][processor][i];
	}

	for ( i = 0; i < num_new_octs; i++ ) {
		new_buffer[ num_local_buffers[level][processor] + i ] = local_octs[i];
	}

	num_local_buffers[level][processor] += num_new_octs;
	cart_free( local_buffers[level][processor] );
	local_buffers[level][processor] = new_buffer;

	cart_assert( level == max_level || num_local_buffers[level][processor] >= num_local_buffers[level+1][processor] / num_children );
}

void cell_buffer_delete_local_oct( int local_oct, int processor ) {
	int i, j;
	int *new_buffer;
	int level;
	int remote_oct;
	int found = 0;

	level = oct_level[local_oct];

	for ( i = 0; i < num_local_buffers[level][processor]; i++ ) {
		if ( local_buffers[level][processor][i] == local_oct ) {
			for ( j = i+1; j < num_local_buffers[level][processor]; j++ ) {
				local_buffers[level][processor][j-1] = local_buffers[level][processor][j];
			}
			found = 1;
			break;
		}
	}

	cart_assert( found );

	num_local_buffers[level][processor]--;

	remote_oct = index_hash_lookup( buffer_oct_reverse_hash[processor], local_oct );

	index_hash_delete( buffer_oct_reverse_hash[processor], remote_oct );
	index_hash_delete( buffer_oct_hash[processor], local_oct );

	new_buffer = cart_alloc(int, num_local_buffers[level][processor] );
	for ( i = 0; i < num_local_buffers[level][processor]; i++ ) {
		new_buffer[i] = local_buffers[level][processor][i];
	}
	cart_free( local_buffers[level][processor] );

	local_buffers[level][processor] = new_buffer;
}

void cell_buffer_delete_remote_oct( int level, int processor, int local_oct ) {
	int i, j;
	int *new_buffer;
	int found = 0;

	cart_assert( level > min_level && level <= max_level );
	cart_assert( processor >= 0 && processor < num_procs );
	cart_assert( local_oct != NULL_OCT );

	for ( i = 0; i < num_remote_buffers[level][processor]; i++ ) {
		if ( remote_buffers[level][processor][i] == local_oct ) {
			for ( j = i + 1; j < num_remote_buffers[level][processor]; j++ ) {
				remote_buffers[level][processor][j-1] = remote_buffers[level][processor][j];
			}
			found = 1;
			break;
		}
	}

	if ( !found ) {
		cart_error("cell_buffer_delete_remote_oct: %u was not found in list %u long", local_oct, 
				num_remote_buffers[level][processor] );
	} else {
		num_remote_buffers[level][processor]--;
		new_buffer = cart_alloc(int, num_remote_buffers[level][processor] );
		for ( i = 0; i < num_remote_buffers[level][processor]; i++ ) {
			new_buffer[i] = remote_buffers[level][processor][i];
		}
		cart_free( remote_buffers[level][processor] );
		remote_buffers[level][processor] = new_buffer;
	}
}

int cell_can_prune( int cell, int proc ) {
	int i;
	int neighbor;
	int sfc;
	int neighbors[num_neighbors];
	int neighbors2[num_secondary_neighbors];
	int can_prune = 1;

	cart_assert( cell >= 0 && cell < num_cells );
	cart_assert( cell_is_local(cell) );

	/* check direct neighbors */
	cell_all_neighbors( cell, neighbors );

	if ( proc == PROC_ANY ) {
		for ( i = 0; i < num_neighbors; i++ ) {
			cart_assert( neighbors[i] == NULL_OCT || cell_parent_root_sfc(neighbors[i]) >= 0 );
			if ( neighbors[i] == NULL_OCT ||
					root_cell_type( cell_parent_root_sfc(neighbors[i]) ) == CELL_TYPE_BUFFER ) {
				can_prune = 0;
				break;
			}
		}
	
		/* check secondary neighbors */ 
		for ( i = 0; can_prune && i < num_secondary_neighbors; i++ ) {
			neighbors2[i] = cell_neighbor( neighbors[ secondary_neighbors[i][0] ],
					secondary_neighbors[i][1] );
	
			cart_assert( neighbors2[i] == NULL_OCT || cell_parent_root_sfc(neighbors2[i]) >= 0 );
			if ( neighbors2[i] == NULL_OCT ||
					root_cell_type( cell_parent_root_sfc(neighbors2[i]) ) == CELL_TYPE_BUFFER ) {
				can_prune = 0;
				break;
			}
		}

		/* check tertiary neighbors */
		for ( i = 0; can_prune && i < num_tertiary_neighbors; i++ ) {
			neighbor = cell_neighbor( neighbors2[ tertiary_neighbors[i][0] ],
					tertiary_neighbors[i][1] );

			cart_assert( neighbor == NULL_OCT || cell_parent_root_sfc(neighbor) >= 0 );
			if ( neighbor == NULL_OCT ||
					root_cell_type( cell_parent_root_sfc(neighbor) ) == CELL_TYPE_BUFFER ) {
				can_prune = 0;
				break;
			}
		}
	} else {
		/* check primary neighbors */
		for ( i = 0; i < num_neighbors; i++ ) {
			if ( neighbors[i] == NULL_OCT ) {
				can_prune = 0;
				break;
			} else {
				sfc = cell_parent_root_sfc( neighbors[i] );
				cart_assert( sfc >= 0 && sfc < max_sfc_index );
	
				if ( processor_owner(sfc) == proc ||
						root_cell_type(sfc) == CELL_TYPE_BUFFER ) {
					can_prune = 0;
					break;
				}
			}
		}

		/* check secondary neighbors */
		for ( i = 0; can_prune && i < num_secondary_neighbors; i++ ) {
			neighbors2[i] = cell_neighbor( neighbors[ secondary_neighbors[i][0] ],
					secondary_neighbors[i][1] );

			if ( neighbors2[i] == NULL_OCT ) {
				can_prune = 0;
				break;
			} else {
				sfc = cell_parent_root_sfc( neighbors2[i] );
				cart_assert( sfc >= 0 && sfc < max_sfc_index );

				if ( processor_owner(sfc) == proc ||
						root_cell_type(sfc) == CELL_TYPE_BUFFER ) {
					can_prune = 0;
					break;
				}
			}
		}

		/* check tertiary neighbors */
		for ( i = 0; can_prune && i < num_tertiary_neighbors; i++ ) {
			neighbor = cell_neighbor( neighbors2[ tertiary_neighbors[i][0] ],
					tertiary_neighbors[i][1] );

			if ( neighbor == NULL_OCT ) {
				can_prune = 0;
				break;
			} else {
				sfc = cell_parent_root_sfc( neighbor );
				cart_assert( sfc >= 0 && sfc < max_sfc_index );

				if ( processor_owner(sfc) == proc ||
						root_cell_type(sfc) == CELL_TYPE_BUFFER ) {
					can_prune = 0;
					break;
				}
			}
		}
	}

	return can_prune;
}

void build_root_cell_buffer() {
	int i;
	int sfc;
	int neighbors[num_stencil];
	skiplist *buffer_list;
	int *buffer_indices;
	int index;

	if ( num_procs > 1 ) {
		buffer_list = skiplist_init();

		for ( sfc = proc_sfc_index[local_proc_id]; sfc < proc_sfc_index[local_proc_id+1]; sfc++ ) {	
			root_cell_uniform_stencil( sfc, neighbors );

			for ( i = 0; i < num_stencil; i++ ) {
				if ( !root_cell_is_local(neighbors[i]) ) {
					skiplist_insert( buffer_list, neighbors[i] );
				}
			}
		}

		num_buffer_cells[min_level] = skiplist_size( buffer_list );

		cart_debug("num_buffer_cells[min_level] = %d", num_buffer_cells[min_level] );
		cart_debug("num_cells_per_level[min_level] = %d", num_cells_per_level[min_level] );

		if ( num_buffer_cells[min_level] + num_cells_per_level[min_level] > num_cells ) {
			cart_error("Ran out of cells allocating root buffer");
		}

		buffer_cell_sfc_index = cart_alloc(int, num_buffer_cells[min_level] );
		buffer_indices = cart_alloc(int, num_buffer_cells[min_level] );

		index = 0;
		skiplist_iterate( buffer_list );
		while ( skiplist_next( buffer_list, &sfc ) ) {
			buffer_indices[index] = index + num_cells_per_level[min_level];
			buffer_cell_sfc_index[index] = sfc;
			index++;
		}

		cart_assert( index == num_buffer_cells[min_level] );
		skiplist_destroy( buffer_list );

		buffer_root_hash = index_hash_create( num_buffer_cells[min_level] );
		index_hash_add_list( buffer_root_hash, num_buffer_cells[min_level], 
				buffer_cell_sfc_index, buffer_indices );

		cart_free( buffer_indices );
	} else {
		num_buffer_cells[min_level] = 0;
	}

	root_buffer_enabled = 1;
}

/*******************************************************
 * build_cell_buffer
 *******************************************************/
void build_cell_buffer()
/* purpose: builds a list of the cells the local processor needs
 *      to have a complete buffer, then all processors communicate
 *      so each processor has a local copy of the needed cells.
 */
{
	int sfc, i;
	int proc, level;
	int processor;
	int neighbors[num_stencil];
	int count[max_level-min_level+1];
	int ioct;
	int num_hash_octs;

	pack *buffer_list;

	cart_assert( root_buffer_enabled == 1 );

	if ( num_procs == 1 ) {
		buffer_enabled = 1;
		return;
	}

	start_time( BUILD_CELL_BUFFER_TIMER );

	cart_assert( buffer_enabled == 0 );

	cart_debug("building cell buffer");

	/* create pack to hold buffer trees */
	buffer_list = pack_init( CELL_TYPE_BUFFER );

	/* loop over all of our local cells, and determine if another processor
	 * needs to buffer them, assumes the stencil is symmetric */
	for ( sfc = proc_sfc_index[local_proc_id]; sfc < proc_sfc_index[local_proc_id+1]; sfc++ ) {
		cart_assert( root_cell_is_local( sfc ) );

		root_cell_uniform_stencil( sfc, neighbors );

		for ( i = 0; i < num_stencil; i++ ) {
			if ( !root_cell_is_local( neighbors[i] ) ) {
				processor = processor_owner( neighbors[i] );
				pack_add_root_tree( buffer_list, processor, sfc );
			}
		}
	}

	pack_apply( buffer_list );

	for ( level = min_level+1; level <= max_level; level++ ) {
		num_buffer_cells[level] = 0;
	}

	/* process received buffer cells */
	for ( proc = 0; proc < num_procs; proc++ ) {
		num_local_buffers[min_level][proc] = buffer_list->num_receiving_cells[proc][min_level];
		num_remote_buffers[min_level][proc] = buffer_list->num_sending_cells[proc][min_level];

		for ( level = min_level+1; level <= max_level; level++ ) {
			num_remote_buffers[level][proc] = buffer_list->num_sending_cells[proc][level] / num_children;
			num_local_buffers[level][proc] = buffer_list->num_receiving_cells[proc][level] / num_children;
		}

		for ( level = min_level; level <= max_level; level++ ) {
			cart_assert( num_remote_buffers[level][proc] >= 0 );
			remote_buffers[level][proc] = cart_alloc(int, num_remote_buffers[level][proc] );
			count[level] = 0;

			local_buffers[level][proc] = cart_alloc(int, num_local_buffers[level][proc] );
		}

		for ( i = 0; i < buffer_list->num_sending_cells[proc][min_level]; i++ ) {
			remote_buffers[min_level][proc][i] = buffer_list->root_cells[proc][i];
		}

		/* pack remote octs */
		for ( i = 0; i < buffer_list->num_sending_cells_total[proc]; i++ ) {
			ioct = buffer_list->cell_refined[proc][i];
			cart_assert( ioct == NULL_OCT || ( ioct >= 0 && ioct < num_octs ) );

			if ( ioct != NULL_OCT ) {
				level = oct_level[ioct];
				remote_buffers[level][proc][count[level]] = ioct;
				count[level]++;
			}			
		}

		/* create hash tables to hold octs */
		num_hash_octs = ( buffer_list->num_sending_cells_total[proc] - 
				buffer_list->num_sending_cells[proc][min_level] ) / num_children;

		buffer_oct_hash[proc] = index_hash_create( num_hash_octs );
		buffer_oct_reverse_hash[proc] = index_hash_create( num_hash_octs );
	}

	pack_communicate( buffer_list );

	pack_destroy( buffer_list );

	buffer_enabled = 1;

	end_time( BUILD_CELL_BUFFER_TIMER );

	cart_debug("buffer enabled");
}

/*******************************************************
 * destroy_cell_buffer
 *******************************************************/
void destroy_cell_buffer() 
/* purpose: deallocates all root cell space taken up by
 * 	cell buffer and update lists 
 */ 
{
	int num_level_cells;
	int *level_cells;
	int i, j;

	if ( num_procs == 1 && num_buffer_cells[min_level] == 0 ) {
		buffer_enabled = 0;
		root_buffer_enabled = 0;
		return;
	}

	/* free all buffer trees */
	select_level( min_level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		cell_free( level_cells[i] );
	}
	cart_free( level_cells );

	cart_free( buffer_cell_sfc_index );
	num_buffer_cells[min_level] = 0;
	buffer_enabled = 0;
	root_buffer_enabled = 0;

	/* clear linked lists */
	for ( i = min_level; i <= max_level; i++ ) {	
		num_buffer_cells[i] = 0;
		buffer_oct_list[i] = NULL_OCT;

		for ( j = 0; j < num_procs; j++ ) {
			cart_free( remote_buffers[i][j] );
			num_remote_buffers[i][j] = 0;
			num_local_buffers[i][j] = 0;
			cart_free( local_buffers[i][j] );
		}
	}

	for ( i = 0; i < num_procs; i++ ) {
		if ( buffer_oct_hash[i] != NULL ) {
			index_hash_free( buffer_oct_hash[i] );
		}

		if ( buffer_oct_reverse_hash[i] != NULL ) {
			index_hash_free( buffer_oct_reverse_hash[i] );
		}
	}

	index_hash_free( buffer_root_hash );

	cart_debug("cell buffer disabled");
}

/*******************************************************
 * update_buffer_level
 *******************************************************/
void update_buffer_level( int level, const int *var_indices, int num_update_vars ) 
/* purpose: updates all remotely buffered cells of the given
 * 	level and receives updates for all locally buffered
 * 	cells
 */
{
	int i, j, k;
	int var_counter;
	int index;
	int icell;
	int proc;
	int buffer_size, buffer_ptr, recv_count;

	MPI_Request requests[MAX_PROCS];
	MPI_Request receives[MAX_PROCS];
	MPI_Status statuses[MAX_PROCS];
	MPI_Status status;

	float *buffer;
	int recv_offset[MAX_PROCS];

	start_time( COMMUNICATION_TIMER );
	start_time( UPDATE_TIMER );

	cart_assert( num_update_vars > 0 && num_update_vars <= num_vars );
	cart_assert( level >= min_level && level <= max_level );

	/* allocate buffer */
	buffer_size = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		cart_assert( num_local_buffers[level][proc] >= 0 );
		cart_assert( num_remote_buffers[level][proc] >= 0 );

		if ( level == min_level ) {
			buffer_size += num_update_vars * 
				( num_local_buffers[min_level][proc] + num_remote_buffers[min_level][proc] );
		} else {
			buffer_size += num_update_vars * num_children *
				( num_local_buffers[level][proc] + num_remote_buffers[level][proc] );	
		}
	}

	buffer = cart_alloc(float, buffer_size );

	/* set up receives */
	buffer_ptr = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_local_buffers[level][proc] > 0 ) {
			recv_offset[proc] = buffer_ptr;

			if ( level == min_level ) {
				recv_count = num_update_vars * num_local_buffers[min_level][proc];
			} else {
				recv_count = num_update_vars * num_children * num_local_buffers[level][proc];
			}

			MPI_Irecv( &buffer[buffer_ptr], recv_count, MPI_FLOAT, proc, num_update_vars,
				MPI_COMM_WORLD, &receives[proc] );

			buffer_ptr += recv_count;
		} else {
			receives[proc] = MPI_REQUEST_NULL;
		}
	}

	/* set up sends */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_remote_buffers[level][proc] > 0 ) {
			var_counter = buffer_ptr;

			if ( level == min_level ) {
				for ( i = 0; i < num_remote_buffers[min_level][proc]; i++ ) {
					cart_assert( remote_buffers[min_level][proc][i] >= 0 && 
							remote_buffers[min_level][proc][i] < max_sfc_index );

					index = root_cell_location( remote_buffers[min_level][proc][i] );
					for ( j = 0; j < num_update_vars; j++ ) {
						buffer[var_counter++] = cell_var( index, var_indices[j] );
					}
				}

				MPI_Isend( &buffer[buffer_ptr], num_remote_buffers[min_level][proc] * num_update_vars, 
					MPI_FLOAT, proc, num_update_vars, MPI_COMM_WORLD, &requests[proc] );
			} else {
				for ( i = 0; i < num_remote_buffers[level][proc]; i++ ) {
					for ( j = 0; j < num_children; j++ ) {
						icell = oct_child( remote_buffers[level][proc][i], j );
						for ( k = 0; k < num_update_vars; k++ ) {
							buffer[var_counter++] = cell_var( icell, var_indices[k] );
						}
					}
				}
			
				MPI_Isend( &buffer[buffer_ptr], num_remote_buffers[level][proc] * num_update_vars * num_children, 
					MPI_FLOAT, proc, num_update_vars, MPI_COMM_WORLD, &requests[proc] );
			}

			buffer_ptr = var_counter;
		} else {
			requests[proc] = MPI_REQUEST_NULL;
		}
	}

	cart_assert( buffer_size == buffer_ptr );

	start_time( UPDATE_RECV_TIMER );

	do {
		MPI_Waitany( num_procs, receives, &proc, &status );

		if ( proc != MPI_UNDEFINED ) {
			var_counter = recv_offset[proc];
			if ( level == min_level ) {
				for ( i = 0; i < num_local_buffers[min_level][proc]; i++ ) {
					index = local_buffers[min_level][proc][i];

					for ( j = 0; j < num_update_vars; j++ ) {
						cell_var( index, var_indices[j] ) = buffer[var_counter++];
					}
				}
			} else {
				for ( i = 0; i < num_local_buffers[level][proc]; i++ ) {
					index = local_buffers[level][proc][i];

					cart_assert( index >= 0 && index < num_octs );
					cart_assert( !cell_is_local( oct_parent_cell[index] ) );
					cart_assert( processor_owner( oct_parent_root_sfc[index] ) == proc );
					
					for ( j = 0; j < num_children; j++ ) {
						icell = oct_child( index, j );
						for ( k = 0; k < num_update_vars; k++ ) {
							cell_var( icell, var_indices[k] ) = buffer[var_counter++];
						}
					}
				}
			}
		}
	} while ( proc != MPI_UNDEFINED );

	end_time( UPDATE_RECV_TIMER );

	/* wait for all sends to complete (later add latency hiding, delay freeing buffers?) */
	start_time( UPDATE_SEND_TIMER );
	MPI_Waitall( num_procs, requests, statuses );
	end_time( UPDATE_SEND_TIMER );

	cart_free(buffer);

	end_time( UPDATE_TIMER );
	end_time( COMMUNICATION_TIMER );
}

#ifdef PARTICLES
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)

void merge_buffer_cell_densities( int level ) {
	int i;
	int index, child;
	int icell, proc;

#ifdef RT_VAR_SOURCE
	const int vars_per_cell = 3;
	const int density_vars_size = 2;
	const int density_vars[2] = { VAR_DENSITY, RT_VAR_SOURCE };
#else
	const int vars_per_cell = 2;
	const int density_vars_size = 1;
	const int density_vars[1] = { VAR_DENSITY };
#endif

	MPI_Request sends[MAX_PROCS];
	MPI_Request receives[MAX_PROCS];
	MPI_Status status;
	MPI_Status statuses[MAX_PROCS];

	float *send_buffer;
	float *recv_buffer;
	int recv_offset[MAX_PROCS];
	int num_send_vars[MAX_PROCS];
	int num_recv_vars[MAX_PROCS];
	int total_send_vars, total_recv_vars;
	int send_offset;
	int send_count, recv_count;

	start_time( MERGE_DENSITY_TIMER );
	start_time( COMMUNICATION_TIMER );

	cart_assert( buffer_enabled );

	total_send_vars = total_recv_vars = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		cart_assert( num_remote_buffers[level][proc] >= 0 );
		cart_assert( num_local_buffers[level][proc] >= 0 );

		if ( level == min_level ) {
			num_recv_vars[proc] = vars_per_cell*num_remote_buffers[min_level][proc];
			num_send_vars[proc] = vars_per_cell*num_local_buffers[min_level][proc];
		} else {
			num_recv_vars[proc] = vars_per_cell*num_children*num_remote_buffers[level][proc];
			num_send_vars[proc] = vars_per_cell*num_children*num_local_buffers[level][proc];
		}

		cart_assert( num_send_vars[proc] >= 0 && num_recv_vars[proc] >= 0 );

		total_send_vars += num_send_vars[proc];
		total_recv_vars += num_recv_vars[proc];
	}

	/* set up receives */
	recv_buffer = cart_alloc(float, total_recv_vars );

	recv_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_recv_vars[proc] > 0 ) {
			recv_offset[proc] = recv_count;
			MPI_Irecv( &recv_buffer[recv_count], num_recv_vars[proc], MPI_FLOAT,
				proc, 0, MPI_COMM_WORLD, &receives[proc] );
			recv_count += num_recv_vars[proc];
		} else {
			receives[proc] = MPI_REQUEST_NULL;
		}
	}	

	/* pack cell ids and densities */
	send_buffer = cart_alloc(float, total_send_vars );

	send_offset = send_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_send_vars[proc] > 0 ) {
			if ( level == min_level ) {
				cart_assert( num_send_vars[proc] == vars_per_cell*num_local_buffers[level][proc] );
				for ( i = 0; i < num_local_buffers[min_level][proc]; i++ ) {
					icell = local_buffers[min_level][proc][i];

					send_buffer[send_count++] = cell_density(icell);
					send_buffer[send_count++] = cell_first_species_mass(icell); 
#ifdef RT_VAR_SOURCE
					send_buffer[send_count++] = cell_rt_source(icell);
#endif
				}
			} else {
				cart_assert( num_send_vars[proc] == vars_per_cell*num_children*num_local_buffers[level][proc] );

				for ( i = 0; i < num_local_buffers[level][proc]; i++ ) {
					index = local_buffers[level][proc][i];
					cart_assert( index >= 0 && index < num_octs );
					cart_assert( oct_level[index] == level );

					for ( child = 0; child < num_children; child++ ) {
						icell = oct_child( index, child );

						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						send_buffer[send_count++] = cell_density(icell);
						send_buffer[send_count++] = cell_first_species_mass(icell);
#ifdef RT_VAR_SOURCE
						send_buffer[send_count++] = cell_rt_source(icell);
#endif
					}
				}
			}

#ifdef DEBUG
			if ( send_offset + num_send_vars[proc] != send_count ) {
				for ( proc = 0; proc < num_procs; proc++ ) {
					cart_debug("proc = %d, num_send = %d, num_local = %d", 
						proc, num_send_vars[proc], num_local_buffers[level][proc] );
				}
				cart_debug("level = %d", level );
				cart_debug("total_send_vars = %d", total_send_vars );
				cart_debug("i = %d", i );
				cart_debug("send_offset = %d", send_offset );
				cart_debug("num_send_vars[%u] = %d", proc, num_send_vars[proc] );
				cart_debug("num_local_buffers[%d][%d] = %d", level, proc );
				cart_debug("send_count = %d", send_count );
			}
#endif

			cart_assert( send_count <= total_send_vars );
			cart_assert( send_offset + num_send_vars[proc] == send_count );

			MPI_Isend( &send_buffer[send_offset], num_send_vars[proc], MPI_FLOAT,
				proc, 0, MPI_COMM_WORLD, &sends[proc] );

			send_offset = send_count;
		} else {
			sends[proc] = MPI_REQUEST_NULL;
		}
	}

	/* process receives as they come in */
	do {
		MPI_Waitany( num_procs, receives, &proc, &status );

		if ( proc != MPI_UNDEFINED ) {
			recv_count = recv_offset[proc];

			if ( level == min_level ) {
				for ( i = 0; i < num_remote_buffers[min_level][proc]; i++ ) {
					icell = root_cell_location( remote_buffers[min_level][proc][i] );

					cell_density(icell) += recv_buffer[recv_count++];
					cell_first_species_mass(icell) += recv_buffer[recv_count++];
#ifdef RT_VAR_SOURCE
					cell_rt_source(icell) += recv_buffer[recv_count++];
#endif
				}
			} else {
				for ( i = 0; i < num_remote_buffers[level][proc]; i++ ) {
					index = remote_buffers[level][proc][i];

					for ( child = 0; child < num_children; child++ ) {
						icell = oct_child( index, child );

						cell_density(icell) += recv_buffer[recv_count++];
						cell_first_species_mass(icell) += recv_buffer[recv_count++];
#ifdef RT_VAR_SOURCE
						cell_rt_source(icell) += recv_buffer[recv_count++];
#endif
					}
				}
			}

			cart_assert( recv_offset[proc] + num_recv_vars[proc] == recv_count );
		}
	} while ( proc != MPI_UNDEFINED );

	cart_free( recv_buffer );

	end_time(COMMUNICATION_TIMER );

	/* now update density variables */
	start_time( MERGE_DENSITIES_UPDATE_TIMER );
	update_buffer_level( level, density_vars, density_vars_size );
	end_time( MERGE_DENSITIES_UPDATE_TIMER );

	/* wait for sends */
	start_time( COMMUNICATION_TIMER );
	MPI_Waitall( num_procs, sends, statuses );
	end_time( COMMUNICATION_TIMER );

	cart_free( send_buffer );

	end_time( MERGE_DENSITY_TIMER );
}

#endif /* GRAVITY */
#endif /* PARTICLES */

void split_buffer_cells( int level, int *cells_to_split, int num_cells_to_split ) {
	int i, j, k, m;
	int proc;
	int buffer_size;
	int info_per_cell;
	int root, child, ichild, parent;
	int cell_index, oct_index;
	int result;
	int *buffer_cells_to_split[MAX_PROCS];
	float *new_buffer_cell_vars;
	int *cells_split[MAX_PROCS];
	float *new_cell_vars[MAX_PROCS];
	int *remote_buffer_list;
	int *remote_octs;
	int *local_octs;
	int *new_octs;
	int num_new_octs;
	int num_cells_to_send;
	int num_cells_to_recv;
	int num_cells_split;
	int num_cell_vars;
	int new_cell_count;

	MPI_Status status;
	MPI_Request receives[MAX_PROCS];
	MPI_Request sends[2*MAX_PROCS];

	cart_assert( level < max_level );

	start_time( COMMUNICATION_TIMER );

	/* need 2 ints to describe split root cell, 3 for all other levels */
	if ( level == min_level ) {
		info_per_cell = 2;
	} else {
		info_per_cell = 3;
	}

	/* set up receive buffers */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_local_buffers[level][proc] > 0 ) {
			/* cells already split cannot be split, so subtract # of level+1 octs */
			if ( level == min_level ) {
				buffer_size = num_local_buffers[level][proc] - num_local_buffers[level+1][proc];
			} else {
				buffer_size = num_children*num_local_buffers[level][proc] - num_local_buffers[level+1][proc];
			}

			cart_assert( buffer_size >= 0 );

			buffer_cells_to_split[proc] = cart_alloc(int, info_per_cell*buffer_size );

			MPI_Irecv( buffer_cells_to_split[proc], info_per_cell*buffer_size, MPI_INT,
					proc, 0, MPI_COMM_WORLD, &receives[proc] );
		} else {
			receives[proc] = MPI_REQUEST_NULL;
		}			
	}

	/* generate list of cells to send to each proc */

	/* prune split cells to those which may be buffered */
	for ( i = 0; i < num_cells_to_split; i++ ) {
		if ( cells_to_split[i] != NULL_OCT && cell_can_prune( cells_to_split[i], PROC_ANY ) ) {
			cells_to_split[i] = -1;
		}
	}

	/* move -1's to the end */
	qsort( cells_to_split, num_cells_to_split, sizeof(int), compare_ints );
	
	/* now how many do we have left? */
	for ( num_cells_split = 0; num_cells_split < num_cells_to_split; num_cells_split++ ) {
		if ( cells_to_split[num_cells_split] == -1 ) {
			break;
		}
	}

	/* count number split by each proc */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_remote_buffers[level][proc] > 0 ) {
			num_cells_to_send = 0;

			/* create and sort remote_buffer_list (since order in remote_buffers must be preserved) */
			remote_buffer_list = cart_alloc(int, num_remote_buffers[level][proc] );

			for ( i = 0; i < num_remote_buffers[level][proc]; i++ ) {
				remote_buffer_list[i] = remote_buffers[level][proc][i];
			}

			qsort( remote_buffer_list, num_remote_buffers[level][proc], sizeof(int), compare_ints );

			/* loop through split cells and count those buffered by proc
			 * assumes remote_buffers and cells_to_split are both sorted */
			k = 0;
			if ( level == min_level ) {
				for ( j = 0; j < num_cells_split; j++ ) {
					root = cell_parent_root_sfc( cells_to_split[j] );
					while ( (k < num_remote_buffers[min_level][proc]) && 
							(root > remote_buffer_list[k]) ) {
						k++;
					}
					
					if ( k < num_remote_buffers[min_level][proc]
							&& root == remote_buffer_list[k]
							&& !cell_can_prune( cells_to_split[j], proc ) ) {

						num_cells_to_send++;
					}
				}
			} else {
				for ( j = 0; j < num_cells_split; j++ ) {
					parent = cell_parent_oct( cells_to_split[j] );
					while ( k < num_remote_buffers[level][proc] && 
							parent > remote_buffer_list[k] ) {
						k++;
					}

					if ( k < num_remote_buffers[level][proc]
							&& parent == remote_buffer_list[k]
							&& !cell_can_prune( cells_to_split[j], proc ) ) {

						num_cells_to_send++;
					}
				}
			}

			/* now allocate space to store these cells */
			cells_split[proc] = cart_alloc(int, info_per_cell*num_cells_to_send );
			new_cell_vars[proc] = cart_alloc(float, num_vars*num_children*num_cells_to_send );

			/* save up oct list to add to remote_buffers array all at once */
			new_octs = cart_alloc(int, num_cells_to_send );
			num_new_octs = 0;

			new_cell_count = 0;
			num_cell_vars = 0;

			/* loop again through split cells and pack the new cell variables */
			k = 0;
			if ( level == min_level ) {
				for ( j = 0; j < num_cells_split; j++ ) {
					root = cell_parent_root_sfc( cells_to_split[j] );
					while ( (k < num_remote_buffers[min_level][proc]) && 
							(root > remote_buffer_list[k]) ) {
						k++;
					}

					if ( k < num_remote_buffers[min_level][proc]
							&& root == remote_buffer_list[k]
							&& !cell_can_prune( cells_to_split[j], proc ) ) {

						cells_split[proc][new_cell_count++] = root;
						cells_split[proc][new_cell_count++] = cell_child_oct[ cells_to_split[j] ];

						new_octs[num_new_octs++] = cell_child_oct[ cells_to_split[j] ];

						/* pack each new child */
						for ( child = 0; child < num_children; child++ ) {
							ichild = cell_child( cells_to_split[j], child );
							for ( m = 0; m < num_vars; m++ ) {
								new_cell_vars[proc][num_cell_vars++] = cell_var(ichild,m);
							}
						}
					}
				}
			} else {
				for ( j = 0; j < num_cells_split; j++ ) {
					parent = cell_parent_oct( cells_to_split[j] );
					while ( k < num_remote_buffers[level][proc] && parent > remote_buffer_list[k] ) {
						k++;
					}

					if ( k < num_remote_buffers[level][proc]
							&& parent == remote_buffer_list[k]
							&& !cell_can_prune( cells_to_split[j], proc ) ) {

						cells_split[proc][new_cell_count++] = parent;
						cells_split[proc][new_cell_count++] = cell_child_number( cells_to_split[j] );
						cells_split[proc][new_cell_count++] = cell_child_oct[ cells_to_split[j] ];

						new_octs[num_new_octs++] = cell_child_oct[ cells_to_split[j] ];

						for ( child = 0; child < num_children; child++ ) {
							ichild = cell_child( cells_to_split[j], child );
							for ( m = 0; m < num_vars; m++ ) {
								new_cell_vars[proc][num_cell_vars++] = cell_var(ichild,m);
							}
						}
					}
				}
			}

			cart_assert( new_cell_count == info_per_cell * num_cells_to_send );
			cart_assert( num_cell_vars == num_vars * num_children * num_cells_to_send );

			/* send split cell arrays to other proc */
			MPI_Isend( cells_split[proc], info_per_cell*num_cells_to_send, MPI_INT, 
				proc, 0, MPI_COMM_WORLD, &sends[proc] );
			MPI_Isend( new_cell_vars[proc], num_cell_vars, MPI_FLOAT, proc, 0, 
				MPI_COMM_WORLD, &sends[num_procs+proc] );

			/* add all newly buffered octs to remote_buffers */
			cell_buffer_add_remote_octs( level+1, proc, new_octs, num_new_octs );
			cart_free( new_octs );

			cart_free( remote_buffer_list );
                } else {
                        sends[proc] = MPI_REQUEST_NULL;
                        sends[num_procs+proc] = MPI_REQUEST_NULL;
                }
	}

	/* now wait to receive from each proc */
	do {
		MPI_Waitany( num_procs, receives, &proc, &status );

                if ( proc != MPI_UNDEFINED ) {
                        MPI_Get_count( &status, MPI_INT, &buffer_size );

			cart_assert( buffer_size % info_per_cell == 0 );

			num_cells_to_recv = buffer_size / info_per_cell;

			new_buffer_cell_vars = cart_alloc(float, num_children*num_vars*num_cells_to_recv );

			MPI_Recv( new_buffer_cell_vars, num_children*num_vars*num_cells_to_recv, MPI_FLOAT,
				proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

			num_new_octs = 0;
			num_cell_vars = 0;

			remote_octs = cart_alloc(int, num_cells_to_recv );
			local_octs = cart_alloc(int, num_cells_to_recv );

			if ( level == min_level ) {
				for ( j = 0; j < num_cells_to_recv; j++ ) {
					root = root_cell_location( buffer_cells_to_split[proc][2*j] );
					cart_assert( cell_level(root) == min_level );
					cart_assert( !cell_is_local(root) );

					/* split buffer cell */
					result = split_cell( root );
					cart_assert( result == 0 );

					/* copy over cell vars */
					for ( child = 0; child < num_children; child++ ) {
						ichild = cell_child(root,child);
						for ( m = 0; m < num_vars; m++ ) {
							cell_var(ichild,m) = new_buffer_cell_vars[num_cell_vars++];
						}
					}

					remote_octs[num_new_octs] = buffer_cells_to_split[proc][2*j+1];
					local_octs[num_new_octs] = cell_child_oct[root];
					num_new_octs++;
				}
			} else {
				for ( j = 0; j < num_cells_to_recv; j++ ) {
                                        cart_assert( buffer_cells_to_split[proc][3*j] != NULL_OCT );

					oct_index = index_hash_lookup( buffer_oct_hash[proc],
									 buffer_cells_to_split[proc][3*j] );
					cell_index = oct_child( oct_index, buffer_cells_to_split[proc][3*j+1] );

					cart_assert( !cell_is_local(cell_index) );
					cart_assert( cell_level(cell_index) == level );

					result = split_cell( cell_index );
					cart_assert( result == 0 );

					/* copy over cells */
					for ( child = 0; child < num_children; child++ ) {
						ichild = cell_child(cell_index,child);
						for ( m = 0; m < num_vars; m++ ) {
							cell_var(ichild,m) = new_buffer_cell_vars[num_cell_vars++];
						}
					}

					remote_octs[num_new_octs] = buffer_cells_to_split[proc][3*j+2];
					local_octs[num_new_octs] = cell_child_oct[cell_index];
					num_new_octs++;
				}
			}

			index_hash_add_list( buffer_oct_hash[proc], num_new_octs, remote_octs, local_octs );
			index_hash_add_list( buffer_oct_reverse_hash[proc], num_new_octs, local_octs, remote_octs );

			cell_buffer_add_local_octs( level+1, proc, local_octs, num_new_octs );

			cart_free( remote_octs );
			cart_free( local_octs );

			cart_free( buffer_cells_to_split[proc] );
			cart_free( new_buffer_cell_vars );
		}
	} while ( proc != MPI_UNDEFINED );

	/* wait for all sends to complete */
	MPI_Waitall( 2*num_procs, sends, MPI_STATUSES_IGNORE );

	/* free all send buffers */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_remote_buffers[level][proc] > 0 ) {
			cart_free( cells_split[proc] );
			cart_free( new_cell_vars[proc] );
		}
	}

        end_time( COMMUNICATION_TIMER );
}

int *oct_list;
int compare_octs_to_join( const void *a, const void *b ) {
	return ( oct_list[*(int *)a] - oct_list[*(int *)b] );
}

void join_buffer_cells( int level, int *octs_to_join, int *parent_root_sfc, int num_octs_to_join ) {
	int i, j, k;
	int result;
	int *buffer_octs_to_join[MAX_PROCS];
	int *octs_joined[MAX_PROCS];
	int *remote_buffer_list;
	int num_octs_joined;
	int num_octs_to_send;
	int oct_index, cell_index;
	int *order;
	int index;
	MPI_Status statuses[MAX_PROCS];
	MPI_Request receives[MAX_PROCS];
	MPI_Request sends[MAX_PROCS];

	start_time( COMMUNICATION_TIMER );
	
	/* set up receives */
	for ( i = 0; i < num_procs; i++ ) {
		if ( num_local_buffers[level+1][i] > 0 ) {
		  buffer_octs_to_join[i] = cart_alloc(int, num_local_buffers[level+1][i] );
			MPI_Irecv( buffer_octs_to_join[i], num_local_buffers[level+1][i], 
					MPI_INT, i, 0, MPI_COMM_WORLD, &receives[i] );
		} else {
			receives[i] = MPI_REQUEST_NULL;
		}
	}

	order = cart_alloc(int, num_octs_to_join );
	for ( i = 0; i < num_octs_to_join; i++ ) {
		order[i] = i;
	}

	oct_list = octs_to_join;
	qsort( order, num_octs_to_join, sizeof(int), compare_octs_to_join );

	for ( i = 0; i < num_procs; i++ ) {
		if ( num_remote_buffers[level+1][i] > 0 ) {
			octs_joined[i] = cart_alloc(int, num_remote_buffers[level+1][i] );
			num_octs_to_send = 0;

			remote_buffer_list = cart_alloc(int, num_remote_buffers[level+1][i] );
			for ( j = 0; j < num_remote_buffers[level+1][i]; j++ ) {
				remote_buffer_list[j] = remote_buffers[level+1][i][j];
			}

			qsort( remote_buffer_list, num_remote_buffers[level+1][i], sizeof(int), compare_ints );

			/* build list of joined octs buffered by this processor  */
			k = 0;
			for ( j = 0; j < num_octs_to_join; j++ ) {
				index = order[j];

				while ( k < num_remote_buffers[level+1][i] && 
						octs_to_join[index] > remote_buffer_list[k] ) {
					k++;
				}
				if ( k < num_remote_buffers[level+1][i] && 
						octs_to_join[index] == remote_buffer_list[k] ) {
					octs_joined[i][num_octs_to_send++] = octs_to_join[index];
				}
			}

			MPI_Isend( octs_joined[i], num_octs_to_send, MPI_INT, i, 0, MPI_COMM_WORLD, &sends[i] );

			cart_free( remote_buffer_list );

			/* delete joined octs from remote_oct list */
			for ( j = 0; j < num_octs_to_send; j++ ) {
				cell_buffer_delete_remote_oct( level+1, i, octs_joined[i][j] );
			}
		} else {
			sends[i] = MPI_REQUEST_NULL;
		}
	}

	cart_free( order );

	/* wait for all receives to complete */
	/* MPI_Waitall( num_procs, receives, statuses ); */

	for ( i = 0; i < num_procs; i++ ) {
		if ( num_local_buffers[level+1][i] > 0 ) {
			MPI_Wait( &receives[i], &statuses[i] );
			MPI_Get_count( &statuses[i], MPI_INT, &num_octs_joined );

			for ( j = 0; j < num_octs_joined; j++ ) {
				cart_assert( buffer_octs_to_join[i][j] >= 0 && buffer_octs_to_join[i][j] < num_octs );

				oct_index = index_hash_lookup( buffer_oct_hash[i], buffer_octs_to_join[i][j] );
				cart_assert( oct_index >= 0 && oct_index < num_octs );
				cart_assert( oct_level[oct_index] == level+1 );

				cell_index = oct_parent_cell[oct_index];
				cell_buffer_delete_local_oct( oct_index, i );

				cart_assert( cell_is_refined(cell_index) );
				result = join_cell( cell_index );
				cart_assert( result == 0 );
			}

			cart_free( buffer_octs_to_join[i] );
		}
	}	

	/* wait for all sends to complete */
	MPI_Waitall( num_procs, sends, statuses );

	for ( i = 0; i < num_procs; i++ ) {
		if ( num_remote_buffers[level+1][i] > 0 ) {
			cart_free( octs_joined[i] );
		}
	}	

	end_time( COMMUNICATION_TIMER );
}
