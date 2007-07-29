#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "parallel.h"
#include "tree.h"
#include "defs.h"
#include "particle.h"
#include "iterators.h"
#include "cell_buffer.h"
#include "auxiliary.h"
#include "sfc.h"
#include "pack.h"
#include "timing.h"
#include "io.h"
#include "density.h"
#include "timestep.h"
#include "cache.h" 
#include "hydro_tracer.h"
#include "tree_linkedlist.h"
#include "load_balance.h"

float cost_per_cell		= 1.0;
float cost_per_particle		= 0.1;
float est_root_buffer_fraction	= 0.15;
float est_buffer_fraction	= 0.45;
int load_balance_frequency	= 0;

int divide_list_recursive( float *global_work, int *global_cells, int num_root_cells_in_division,
		double total_work, long total_cells, int num_procs_in_division, int first_proc,
		int first_cell_index, int *proc_index ) {

	int i, j, k;
	float work_frac_left;
	double current_work, current_work_left, current_work_right;
	long current_cells, current_cells_left, current_cells_right;
	long max_cells_left, max_cells_right;
	int num_procs_left, num_procs_right;
	int left_failure, right_failure;

	cart_assert( num_root_cells_in_division >= num_procs_in_division );

	if ( num_procs_in_division == 1 ) {
		proc_index[first_proc] = first_cell_index;
		proc_index[first_proc+1] = first_cell_index+num_root_cells_in_division;
		return 0;
	}

	num_procs_left = num_procs_in_division/2;
	num_procs_right = num_procs_in_division - num_procs_left;

	max_cells_left = (long)num_procs_left*(long)num_cells;
	max_cells_right = (long)num_procs_right*(long)num_cells;

	cart_assert( num_procs_left >= 1 );
	cart_assert( num_procs_right >= 1 );

	work_frac_left = (float)(num_procs_left)/num_procs_in_division;

	current_cells = 0;
	for ( j = 0; j < num_root_cells_in_division; j++ ) {
		if ( 		j >= num_procs_left && 
				( num_root_cells_in_division - j < num_procs_right ||
				j > (long)num_procs_left*(long)((1.0-est_root_buffer_fraction)*(double)num_particles) || 
				(float)j*(1.+est_root_buffer_fraction)+
				(float)(current_cells-j)*(1.+est_buffer_fraction) > max_cells_left ) ) {
			break;
		}
		current_cells += global_cells[j];
	}

	cart_assert( j >= num_procs_left );

	current_cells = 0;
	for ( i = num_root_cells_in_division-1; i >= 0; i-- ) {
		if ( 		num_root_cells_in_division - i > num_procs_right &&
				( i <= num_procs_left ||
				(num_root_cells_in_division - i - 1) > 
				(long)num_procs_right*(long)((1.0-est_root_buffer_fraction)*(double)num_particles ) || 
				(float)( num_root_cells_in_division - i - 1)*(1.+est_root_buffer_fraction) +
				(float)( current_cells - num_root_cells_in_division + i + 1 )*(1.+est_buffer_fraction) >
				max_cells_right ) ) {
			break;
		}
		current_cells += global_cells[i];
	}

	current_work_left = 0.0;
	current_cells_left = 0;

	for ( k = 0; k <= i; k++ ) {
		current_work_left += global_work[k];
		current_cells_left += global_cells[k];
	}
	
	proc_index[first_proc] = first_cell_index;

	if ( i < j ) {
		/* pick the best guess splitting point */
		current_work = 0.0;
		current_cells = 0;

		for ( k = i+1; k < j - 1; k++ ) {
			if ( current_work + current_work_left + global_work[k] > work_frac_left*total_work ) {
				break;
			} else {
				current_work += global_work[k];
				current_cells += global_cells[k];
			}
		}

		if ( num_procs_left > 1 ) {
			left_failure = divide_list_recursive( global_work, global_cells, 
					k, current_work+current_work_left, 
					current_cells+current_cells_left, num_procs_left, first_proc, 
					first_cell_index, proc_index );
		} else {
			left_failure = 0;
			proc_index[first_proc+num_procs_left] = k+first_cell_index;
		}

		if ( num_procs_right > 1 ) {
			right_failure = divide_list_recursive( &global_work[k], &global_cells[k], 
					num_root_cells_in_division - k,
					total_work - current_work - current_work_left, 
					total_cells - current_cells - current_cells_left, num_procs_right, 
					first_proc+num_procs_left, first_cell_index+k, proc_index );
		} else {
			right_failure = 0;
			proc_index[first_proc+num_procs_in_division] = first_cell_index+num_root_cells_in_division;
		}

		/* if we're unable to divide recursively, try linearly instead (for speed) */
		if ( left_failure && right_failure ) {
			return divide_list_linear( global_work, global_cells, num_root_cells_in_division,
					total_work, total_cells, num_procs_in_division, first_proc,
					first_cell_index, proc_index );
		} else if ( left_failure ) {
			return divide_list_linear( global_work, global_cells, k, current_work+current_work_left,
					current_cells+current_cells_left, num_procs_left, first_proc,
					first_cell_index, proc_index );
		} else if ( right_failure ) {
			return divide_list_linear( &global_work[k], &global_cells[k], num_root_cells_in_division - k,
					total_work - current_work - current_work_left,
					total_cells - current_cells - current_cells_left, num_procs_right,
					first_proc+num_procs_left, first_cell_index+k, proc_index );
		} else {
			return 0;
		}
	}

	cart_debug("failed because i >= j, %d %d", i, j);
	return -1;
}

int divide_list_linear( float *global_work, int *global_cells, int num_root_cells_in_division,
                double total_work, long total_cells, int num_procs_in_division, int first_proc,
                int first_cell_index, int *proc_index ) {

	int i;
	int proc;
	int index;
	int min_assigned;
	double local_work;
	long local_cells;
	double ideal_work_per_proc;
	int num_procs_remaining;
	int num_root_cells_remaining;

	cart_debug("dividing linearly for %u procs, first_proc = %u", num_procs_in_division, first_proc );

	ideal_work_per_proc = total_work / (double)num_procs_in_division;

	index = 0;
	num_procs_remaining = num_procs_in_division-1;
	num_root_cells_remaining = num_root_cells_in_division;

	proc_index[first_proc] = first_cell_index;

	for ( proc = first_proc+1; proc < num_procs_in_division + first_proc; proc++ ) {
		/* ensure minimum of # of root cells to ensure no later processor is assigned too many cells */
		min_assigned = max( 1, num_root_cells_remaining - 
				(long)num_procs_remaining*(long)((float)num_cells/
				((float)(1.+est_root_buffer_fraction))));

		local_work = 0.0;
		local_cells = 0;
		proc_index[proc] = proc_index[proc-1];
		for ( i = 0; i < min_assigned; i++ ) {
			local_work += global_work[index];
			local_cells += global_cells[index];
			num_root_cells_remaining--;
			index++;
		}
		proc_index[proc] += min_assigned;

		/* ensure maximum of reduced_root_cells, leaving at least one root cell
		 * for each processor left to come, and not assigning more cells than processor
		 * can hold (max_cells_allowed), as well as allowing to go above ideal_work_per_proc
		 * if rest of procs can't hold the cells remaining
		 */
		while (         index < num_root_cells_in_division - num_procs_remaining &&
				num_root_cells_remaining > num_procs_remaining &&
				proc_index[proc]-proc_index[proc-1] < 
				(int)((1.0-est_root_buffer_fraction)*(float)num_particles) &&
				(float)(proc_index[proc]-proc_index[proc-1])*(1.+est_root_buffer_fraction) + 
				(float)(local_cells - (proc_index[proc]-proc_index[proc-1]))*
				(1.+est_buffer_fraction) < num_cells &&
				( local_work + 0.5*global_work[index] < ideal_work_per_proc ||
				(float)num_root_cells_remaining*(1.+est_root_buffer_fraction) +
				(float)(total_cells-(local_cells-num_root_cells_remaining))*(1.+est_buffer_fraction) >=
				(long)num_procs_remaining*(long)num_cells ) ) {

			local_work += global_work[index];
			local_cells += global_cells[index];
			index++;
			proc_index[proc]++;
			num_root_cells_remaining--;
		}

		/* continually re-evaluate ideal work given what work is left */
		num_procs_remaining--;
		total_work -= local_work;
		total_cells -= local_cells;

		if ( num_procs_remaining > 0 ) {
			ideal_work_per_proc = total_work / (double)(num_procs_remaining);
		}
	}

	proc_index[first_proc+num_procs_in_division] = first_cell_index + num_root_cells_in_division;

	return 0;
}

void load_balance_entire_volume( float *global_work, int *global_cells, int *new_proc_sfc_index ) {
	int i, j;
	int last, count;
	int num_blocks;
	int flag, ret;
	int current_proc, current_root_index;
	int num_reserved_procs;
	int num_procs_in_block;
	int num_root_cells_in_block;
	long total_cells, current_cells;
	double total_work, current_work;
	double ideal_work_per_proc;

	/* mark large blocks of unrefined cells */
	i = 0;
	count = 0;
	last = -1;
	num_blocks = 0;
	num_reserved_procs = 0;
	while ( i < num_root_cells ) {
		if ( global_cells[i] <= 9 ) {
			count++;
		} else {
			/* see how big the refined region is */
			j = i;
			while ( j < num_root_cells && global_cells[j] > 9 ) {
				j++;
			}

			if ( j - i < 50 ) {
				/* mark small refined region as if it were unrefined */
				count += j-i;
				i = j-1;
			} else if ( count > 0.15*min( num_particles, num_cells ) ) {
				/* mark unrefined region for special load allocation */
				for ( j = last+1; j < i; j++ ) {
					global_cells[j] = -global_cells[j];
				}
			
				num_blocks++;

#ifdef PARTICLES
				/* assume ~1 particle per root cell */
				num_reserved_procs += (count-1)/
					(int)((1.0-est_root_buffer_fraction)*(float)num_particles) + 1;
#else
				num_reserved_procs += (count-1)/(int)((1.0-buffer_fraction)*(float)num_cells) + 1;
#endif /* particles */
				count = 0;
				last = i;
			} else {
				count = 0;
				last = i;
			}
		}

		i++;
	}

	cart_debug("num_reserved_procs = %d", num_reserved_procs );

	/* compute total work (now not including large blocks of 0/1-level cells) */
	total_work = 0.0;
	total_cells = 0;
	for ( i = 0; i < num_root_cells; i++ ) {
		if ( global_cells[i] > 0 ) {
			total_work += global_work[i];
			total_cells += global_cells[i];
		}
	}

	new_proc_sfc_index[0] = 0;
	ideal_work_per_proc = total_work / (double)(num_procs - num_reserved_procs);

	cart_debug("total_work = %e, num_reserved_procs = %u, ideal_work_per_proc = %e", total_work,
		num_reserved_procs, ideal_work_per_proc );

	/* allocate in blocks of highly refined cells */
	if ( num_blocks > 0 && num_reserved_procs < num_procs ) {
		last = 0;
		i = 0;
		flag = ( global_cells[0] > 0 ) ? 1 : 0;
		current_proc = 0;
		num_procs_in_block = 0;
		while ( i < num_root_cells ) {
			current_cells = 0;
			current_work = 0.0;

			if ( flag ) {
				while ( i < num_root_cells && global_cells[i] > 0 ) {
					current_cells += global_cells[i];
					current_work += global_work[i];
					i++;
				}

				num_root_cells_in_block = i - last;
				num_procs_in_block = max( 1, (int)floor(current_work / ideal_work_per_proc) );

				num_reserved_procs += num_procs_in_block;
				total_work -= current_work;
				ideal_work_per_proc = total_work / (double)(num_procs - num_reserved_procs);
			} else {
				while ( i < num_root_cells && global_cells[i] < 0 ) {
					global_cells[i] = abs(global_cells[i]);
					current_cells += global_cells[i];
					current_work += global_work[i];
					i++;
				}

				num_root_cells_in_block = i - last;

#ifdef PARTICLES
				/* assume ~1 particle per root cell */
				num_procs_in_block = (num_root_cells_in_block-1)/
					(int)((1.0-est_root_buffer_fraction)*(float)num_particles) + 1;
#else
				num_procs_in_block = (num_root_cells_in_block-1)/
					(int)((1.0-buffer_fraction)*(float)num_cells) + 1;
#endif /* particles */

			}

			num_procs_in_block = min( num_procs_in_block, num_procs - current_proc );

			cart_debug("assigning %u procs %d root cells in %u type block, current_work = %e, current_cells = %d",
					num_procs_in_block, num_root_cells_in_block, flag, current_work, current_cells );

			if ( i < num_root_cells ) {
				flag = ( global_cells[i] < 0 ) ? 0 : 1;

				if ( current_proc + num_procs_in_block == num_procs ) {
					for ( ; i < num_root_cells; i++ ) {
						current_work += global_work[i];

						global_cells[i] = abs(global_cells[i]);
						current_cells += global_cells[i];
					}

					num_root_cells_in_block = num_root_cells - last;
				}
			} else {
				num_procs_in_block = num_procs - current_proc;
			}

			new_proc_sfc_index[current_proc] = last;

			/* try dividing until we get a workable solution */
			do {
				ret = divide_list_recursive( &global_work[last], &global_cells[last],
						num_root_cells_in_block, current_work, 
						current_cells, num_procs_in_block,
						current_proc, last, new_proc_sfc_index );

				if ( ret == -1 ) {
					cart_debug("trying again with more processors");
					num_procs_in_block++;

					if ( current_proc + num_procs_in_block >= num_procs ) {
						cart_error("Unable to find useable configuration");
					}
				}
			} while ( ret == -1 );

			current_proc += num_procs_in_block;
			last = i;
		}
	} else {
		ret = divide_list_recursive( global_work, global_cells, 
				num_root_cells, total_work, total_cells,
				num_procs, MASTER_NODE, 0, new_proc_sfc_index );

		if ( ret == -1 ) {
			cart_error("Unable to find proper load balancing division");
		}
	}

	total_work = 0.0;
	for ( i = 0; i < num_procs; i++ ) {
		current_work = 0.0;
		current_cells = 0;
		for ( j = new_proc_sfc_index[i]; j < new_proc_sfc_index[i+1]; j++ ) {
			current_work += global_work[j];
			current_cells += global_cells[j];
		}

		total_work += current_work;
		cart_debug("new_proc_sfc_index[%u] = %u, %u, work = %e, %u cells", i, 
				new_proc_sfc_index[i], new_proc_sfc_index[i+1], current_work, current_cells );
	}

	cart_debug("total work = %e", total_work );

	/* sanity tests */
	for ( i = 0; i < num_procs; i++ ) {
		cart_assert( new_proc_sfc_index[i] < new_proc_sfc_index[i+1] );
		cart_assert( new_proc_sfc_index[i+1] <= num_root_cells );
	}
}

void load_balance() {
	int i, j;
	int level;
	int root;
	int num_parts;
	float level_cost;
	int proc;
	int first_oct, old_first_oct;
	int next, new_oct;
	int coords[nDim];
	int receive_counts[MAX_PROCS];
	int old_proc_sfc_index[MAX_PROCS+1];
	int new_proc_sfc_index[MAX_PROCS+1];

#ifdef PARTICLES
	int ipart, ipart_next;
	int num_parts_to_send[MAX_PROCS];
	int particle_list_to_send[MAX_PROCS];
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	int tracer;
	int num_tracers_to_send[MAX_PROCS];
	int tracer_list_to_send[MAX_PROCS];
#endif /* HYDRO_TRACERS */

	float *local_work;
	int *local_counts;
	float *global_work;
	int *global_cells;
	int neighbors[num_neighbors];
	int min_assigned;
	double current_work, total_work, ideal_work_per_proc;
	int num_cells_moved;
	int root_cell_shift;
	int *local_index;
	int icell, ioct;
	int num_level_cells;
	int max_cells_allowed;
	int *level_cells;
	int sfc, count;
	int sfc_neighbor;
	int new_num_local_root_cells;
	pack *transfer_cells;

	cart_debug("in load balance");

	/* don't bother doing anything when only 1 processor exists */
	if ( num_procs == 1 ) {
		if ( !buffer_enabled ) {
			build_cell_buffer();
			repair_neighbors();
			cart_debug("done building cell buffer");
		}
		return;
	}

	start_time( LOAD_BALANCE_TIMER );

	/* allocate space to store the work estimators */
	local_work = cart_alloc( num_cells_per_level[min_level]* sizeof(float) );
	local_counts = cart_alloc( num_cells_per_level[min_level]* sizeof(int) );

	#pragma omp parallel for
	for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
		local_work[i] = 0.0;
		local_counts[i] = 0;
	}

	/* compute work per root cell */
	for ( level = min_level; level <= max_level; level++ ) {
		/* give more weight to higher levels of refinement */
		level_cost = (float)(1<<level); 

		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		#pragma omp parallel for private(icell,root,num_parts,ipart)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			root = cell_parent_root_cell( icell );
			sfc = cell_parent_root_sfc( icell );

			local_work[root] += cost_per_cell * level_cost;
			local_counts[root]++;

#ifdef PARTICLES
			/* count number of particles in cell */
			num_parts = 0;
			ipart = cell_particle_list[icell];
			while ( ipart != NULL_PARTICLE ) {
				num_parts++;
				ipart = particle_list_next[ipart];				
			}

			local_work[root] += cost_per_particle * (float)num_parts * level_cost;
#endif /* PARTICLES */
		}

		cart_free( level_cells );
	}

	if ( local_proc_id == MASTER_NODE ) {
		global_work = cart_alloc( num_root_cells * sizeof(float) );
		global_cells = cart_alloc( num_root_cells * sizeof(int) );

		for ( i = 0; i < num_procs; i++ ) {
			receive_counts[i] = proc_sfc_index[i+1] - proc_sfc_index[i];
		}
	}

	/* gather all work to MASTER_NODE */
	MPI_Gatherv( local_work, num_cells_per_level[min_level], MPI_FLOAT,
		global_work, receive_counts, proc_sfc_index, MPI_FLOAT,
		MASTER_NODE, MPI_COMM_WORLD );
	MPI_Gatherv( local_counts, num_cells_per_level[min_level], MPI_INT,
		global_cells, receive_counts, proc_sfc_index, MPI_INT,
		MASTER_NODE, MPI_COMM_WORLD );

	cart_free( local_work );
	cart_free( local_counts );

	/* MASTER_NODE determines new splitting of SFC curve */
	if ( local_proc_id == MASTER_NODE ) {
		load_balance_entire_volume( global_work, global_cells, new_proc_sfc_index );
	}

	/* tell other processors which cells to expect */
	MPI_Bcast( &new_proc_sfc_index, num_procs+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	/* how many do we have now? */
	new_num_local_root_cells = new_proc_sfc_index[local_proc_id+1] - new_proc_sfc_index[local_proc_id];

	/* make sure we free up root octs if new first_oct < old_first_oct */
	old_first_oct = cell_parent_oct( num_cells_per_level[min_level] +
			num_buffer_cells[min_level] ) + 1;

	/* get the cell buffer out of the way (easier to rebuild it later) */
	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}

	/* set up list of cells to send/recv */
	transfer_cells = pack_init( CELL_TYPE_LOCAL );

	/* pack cells to move */
	for ( i = proc_sfc_index[local_proc_id]; i < proc_sfc_index[local_proc_id+1]; i++ ) {
		/* find where root cell is now */
		for ( j = -1; j < num_procs - 1; j++ ) {
			if ( i < new_proc_sfc_index[j+1] ) {
				break;
			}
		}

		cart_assert( j >= 0 && j < num_procs );

		if ( j != local_proc_id ) {
			pack_add_root_tree( transfer_cells, j, i );
		}
	}

	/* finish packing cells so we can move them */
	pack_apply( transfer_cells );

	/* free trees which were just packed */
	for ( proc = 0; proc < num_procs; proc++ ) {
		for ( i = 0; i < transfer_cells->num_sending_cells[proc][min_level]; i++ ) {
			cell_free( root_cell_location( transfer_cells->root_cells[proc][i] ) );
		}
	}

        /* make transfer official */
	for ( i = 0; i <= num_procs+1; i++ ) {
		old_proc_sfc_index[i] = proc_sfc_index[i];
		proc_sfc_index[i] = new_proc_sfc_index[i];
	}

	num_cells_per_level[min_level] = new_num_local_root_cells;

	/* determine how many buffer cells to expect */
        build_root_cell_buffer();

	/* move octs out of the way */
	first_oct = cell_parent_oct( num_cells_per_level[min_level] + 
					num_buffer_cells[min_level] ) + 1;
	next_free_oct = max( next_free_oct, first_oct );

	/* remove all octs < first_oct from free_oct_list */
	ioct = free_oct_list;
	while ( ioct != NULL_OCT ) {
		if ( ioct < first_oct ) {
			if ( ioct == free_oct_list ) {
				free_oct_list = oct_next[ioct];
			} else {
				oct_next[oct_prev[ioct]] = oct_next[ioct];
			}

			if ( oct_next[ioct] != NULL_OCT ) {
				oct_prev[oct_next[ioct]] = oct_prev[ioct];
			}
		}

		ioct = oct_next[ioct];
	}

	/* add old root octs to free_oct_list */
	for ( ioct = first_oct; ioct < old_first_oct; ioct++ ) {
		cart_assert( oct_level[ioct] == FREE_OCT_LEVEL );
		oct_next[ioct] = free_oct_list;
		oct_prev[ioct] = NULL_OCT;

		if ( free_oct_list != NULL_OCT ) {
			oct_prev[free_oct_list] = ioct;
		}
		free_oct_list = ioct;
	}

	/* move non-free octs from root cell space */
	for ( ioct = 0; ioct < first_oct; ioct++ ) {
		if ( oct_level[ioct] != FREE_OCT_LEVEL ) {
			/* move it out of the way */
			new_oct = oct_alloc();
			oct_move( ioct, new_oct );

			/* don't use oct_free on ioct since it'll be placed in free oct list */
                	linked_list_remove( &local_oct_list[oct_level[ioct]], ioct );
			oct_level[ioct] = FREE_OCT_LEVEL;
        		oct_parent_root_sfc[ioct] = NULL_OCT;
			oct_parent_cell[ioct] = NULL_OCT;
        
			for ( i = 0; i < num_neighbors; i++ ) {
				oct_neighbors[ioct][i] = NULL_OCT;
			}
		}
	}

	/* move root trees to proper positions */
	if ( local_proc_id > 0 ) {
		if ( new_proc_sfc_index[local_proc_id] > old_proc_sfc_index[local_proc_id] &&
				new_proc_sfc_index[local_proc_id] < old_proc_sfc_index[local_proc_id+1] ) {

			/* move down */
			num_cells_moved = min( old_proc_sfc_index[local_proc_id+1], new_proc_sfc_index[local_proc_id+1] ) 
						- new_proc_sfc_index[local_proc_id];
			root_cell_shift = new_proc_sfc_index[local_proc_id] - old_proc_sfc_index[local_proc_id];

			for ( i = 0; i < num_cells_moved; i++ ) {
				cell_move( i + root_cell_shift, i );
			}
		} else if ( old_proc_sfc_index[local_proc_id] > new_proc_sfc_index[local_proc_id] &&
			    old_proc_sfc_index[local_proc_id] < new_proc_sfc_index[local_proc_id+1] ) {

			/* move up */
			num_cells_moved = min( old_proc_sfc_index[local_proc_id+1], new_proc_sfc_index[local_proc_id+1] )
						- old_proc_sfc_index[local_proc_id];
			root_cell_shift = old_proc_sfc_index[local_proc_id] - new_proc_sfc_index[local_proc_id];

			for ( i = num_cells_moved-1; i >= 0; i-- ) {
				cell_move( i, i + root_cell_shift );
			}
		}
	}

	/* need to repair neighbors here for can_prune */
	repair_neighbors();

	/* communicate trees to new owners */
	pack_communicate( transfer_cells );

	/* cleanup transfer list */
	pack_destroy( transfer_cells );

	/* reorder cells for cache efficiency */
	cache_reorder_tree(); 

	/* re-enable cell buffer */
	repair_neighbors();
	build_cell_buffer();
	repair_neighbors();

	/* move particles */
#ifdef PARTICLES
	for ( proc = 0; proc < num_procs; proc++ ) {
		num_parts_to_send[proc] = 0;
		particle_list_to_send[proc] = NULL_PARTICLE;
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL ) {
			/* check which processor this particle now belongs to */
			for ( i = 0; i < nDim; i++ ) {
				coords[i] = (int)(particle_x[ipart][i]);
			}

			proc = processor_owner( sfc_index( coords ) );

			if ( proc != local_proc_id ) {
				/* don't need to call delete_particle here since the 
				 * linked lists are in cells which have just been freed */
				particle_level[ipart] = FREE_PARTICLE_LEVEL;
				particle_list_next[ipart] = particle_list_to_send[proc];
				particle_list_to_send[proc] = ipart;
				num_parts_to_send[proc]++;
			}
		}
	}

	/* send and receive particles which switched processor */
	trade_particle_lists( num_parts_to_send, particle_list_to_send, -1 );

	/* reorder particles for cache efficiency */
	cache_reorder_particles();
#endif /* PARTICLES */

	/* move tracers */
#ifdef HYDRO_TRACERS
	for ( proc = 0; proc < num_procs; proc++ ) {
		num_tracers_to_send[proc] = 0;
		tracer_list_to_send[proc] = NULL_TRACER;
	}

	for ( tracer = 0; tracer < num_tracers; tracer++ ) {
		if ( tracer_id[tracer] != NULL_TRACER ) {
			/* check which processor this particle now belongs to */
			for ( i = 0; i < nDim; i++ ) {
				coords[i] = (int)(tracer_x[tracer][i]);
			}

			proc = processor_owner( sfc_index( coords ) );

			if ( proc != local_proc_id ) {
				/* don't need to call delete_tracer here since the 
				 * linked lists are in cells which have just been freed */
				tracer_list_next[tracer] = tracer_list_to_send[proc];
				tracer_list_to_send[proc] = tracer;
				num_tracers_to_send[proc]++;
			}
		}
	}

	/* send and receive particles which switched processor */
	trade_tracer_lists( num_tracers_to_send, tracer_list_to_send, -1 );

	/* reorder tracers? */
#endif /* HYDRO_TRACERS */

	end_time( LOAD_BALANCE_TIMER );
}
