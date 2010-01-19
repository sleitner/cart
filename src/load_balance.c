#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cache.h" 
#include "cell_buffer.h"
#include "control_parameter.h"
#include "density.h"
#include "hydro_tracer.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "pack.h"
#include "parallel.h"
#include "particle.h"
#include "sfc.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "tree_linkedlist.h"


float cost_per_cell		= 1.0;
float cost_per_particle		= 0.25;
float est_buffer_fraction	= 0.5;
int load_balance_frequency	= 0;


void config_init_load_balance()
{
  control_parameter_add2(control_parameter_int,&load_balance_frequency,"frequency:load-balance","load_balance_frequency","frequency (in global time steps) for balancing the load across the nodes. Zero value disables load balancing altogether.");

  control_parameter_add2(control_parameter_float,&cost_per_cell,"cost-per-cell","cost_per_cell","computational cost per cell. This parameter is used in load balancing.");

  control_parameter_add2(control_parameter_float,&cost_per_particle,"cost-per-particle","cost_per_particle","computational cost per particle. This parameter is used in load balancing.");
}


void config_verify_load_balance()
{
  cart_assert(load_balance_frequency >= 0);

  cart_assert(cost_per_cell > 0.0);

  cart_assert(cost_per_particle > 0.0);
}



int divide_list_recursive( float *global_work, 
		int *constrained_quantities, 
		int *per_proc_constraints, 
		int num_root_cells_in_division, double total_work,
		int num_procs_in_division, int first_proc,
		int first_cell_index, int *proc_index ) {

	int i, j, k, c;
	float work_frac_left;
	double current_work, current_work_left;
	long sum_constraints[num_constraints];
	int num_procs_left, num_procs_right;
	int left_failure, right_failure;

	cart_assert( num_root_cells_in_division >= num_procs_in_division );

	if ( num_procs_in_division == 1 ) {
		proc_index[first_proc] = first_cell_index;
		proc_index[first_proc+1] = first_cell_index+num_root_cells_in_division;
		return 0;
	}

	for ( c = 0; c < num_constraints; c++ ) {
		sum_constraints[c] = 0;
	}

	for ( i = 0; i < num_root_cells_in_division; i++ ) {
		for ( c = 0; c < num_constraints; c++ ) {
			sum_constraints[c] += constrained_quantities[num_constraints*i+c];
		}
	}

	for ( c = 0; c < num_constraints; c++ ) {
		if ( sum_constraints[c] > (long)num_procs_in_division*(long)per_proc_constraints[c] ) {
			cart_debug("sum_constraints[%d] = %d vs %d", c, sum_constraints[c],
				(long)num_procs_in_division*(long)per_proc_constraints[c] );
			cart_debug("first_proc = %d", first_proc );
			cart_debug("num_procs_in_division = %d", num_procs_in_division );
			cart_debug("first_cell_index = %d", first_cell_index );
			cart_debug("num_root_cells_in_division = %d", num_root_cells_in_division );
		}
		cart_assert( sum_constraints[c] <= (long)num_procs_in_division*(long)per_proc_constraints[c] );
	}

	num_procs_left = num_procs_in_division/2;
	num_procs_right = num_procs_in_division - num_procs_left;

	cart_assert( num_procs_left >= 1 );
	cart_assert( num_procs_right >= 1 );

	work_frac_left = (float)(num_procs_left)/num_procs_in_division;

	/* find minimal splitting points based on constraints */
	for ( c = 0; c < num_constraints; c++ ) {
		sum_constraints[c] = 0;
	}

	/* use greedy algorithm to find max/min splitting points */
	j = 0;
	k = 0;
	while ( j < num_root_cells_in_division-num_procs_right && k < num_procs_left ) {
		for ( c = 0; c < num_constraints; c++ ) {
			sum_constraints[c] += constrained_quantities[num_constraints*j+c];

			if ( sum_constraints[c] >= per_proc_constraints[c] ) {
				break;
			}
		}

		if ( c < num_constraints ) {
			for ( c = 0; c < num_constraints; c++ ) {
				sum_constraints[c] = 0;
			}
			k++;
		} else {
			j++;
		}
	}

	for ( c = 0; c < num_constraints; c++ ) {
		sum_constraints[c] = 0;
	}

	i = num_root_cells_in_division - 1;
	k = 0;
	while ( i >= num_procs_left && k < num_procs_right ) {
		for ( c = 0; c < num_constraints; c++ ) {
			sum_constraints[c] += constrained_quantities[num_constraints*i+c];

			if ( sum_constraints[c] >= per_proc_constraints[c] ) {
				break;
			}
		}

		if ( c < num_constraints ) {
			for ( c = 0; c < num_constraints; c++ ) {
				sum_constraints[c] = 0;
			}
			k++;
		} else {
			i--;
		}
	}

	current_work_left = 0.0;

	for ( k = 0; k <= i; k++ ) {
		current_work_left += global_work[k];
	}
	
	proc_index[first_proc] = first_cell_index;

	/* k \in (i,j] (k is first cell to the right of division) */
	if ( i < j ) {
		/* pick the best splitting point using work as an estimate */
		current_work = 0.0;

		for ( k = i+1; k < j; k++ ) {
			if ( current_work + current_work_left + global_work[k] > 
					work_frac_left*total_work ) {
				break;
			} else {
				current_work += global_work[k];
			}
		}

		if ( num_procs_left > 1 ) {
			left_failure = divide_list_recursive( global_work, 
					constrained_quantities,
					per_proc_constraints, 
					k, current_work+current_work_left, 
					num_procs_left, first_proc, 
					first_cell_index, proc_index );
		} else {
			left_failure = 0;
			proc_index[first_proc+num_procs_left] = k+first_cell_index;
		}

		if ( num_procs_right > 1 ) {
			right_failure = divide_list_recursive( &global_work[k], 
					&constrained_quantities[num_constraints*k], 
					per_proc_constraints, 
					num_root_cells_in_division - k,
					total_work - current_work - current_work_left, 
					num_procs_right, first_proc+num_procs_left, 
					first_cell_index+k, proc_index );
		} else {
			right_failure = 0;
			proc_index[first_proc+num_procs_in_division] = 
					first_cell_index+num_root_cells_in_division;
		}

		/* if we're unable to divide recursively, try linearly instead (for speed) */
		if ( left_failure && right_failure ) {
			return divide_list_linear( global_work, constrained_quantities, 
					per_proc_constraints, num_root_cells_in_division,
					total_work, num_procs_in_division, first_proc,
					first_cell_index, proc_index );
		} else if ( left_failure || right_failure ) {

			if ( left_failure ) {
				for ( ; i <= k; i++ ) {
					current_work_left += global_work[k];
				}
			} else {
				j = k;
			}

			current_work = 0.0;

			for ( k = i+1; k < j; k += (j-i)/10 ) {
				cart_debug("trying all k... k = %d", k );

				if ( num_procs_left > 1 ) {
					left_failure = divide_list_recursive( global_work,
							constrained_quantities,
							per_proc_constraints,
							k, current_work+current_work_left,
							num_procs_left, first_proc,
							first_cell_index, proc_index );
				} else {
					left_failure = 0;
					proc_index[first_proc+num_procs_left] = k+first_cell_index;
				}

				if ( num_procs_right > 1 ) {
					right_failure = divide_list_recursive( &global_work[k],
							&constrained_quantities[num_constraints*k],
							per_proc_constraints,
							num_root_cells_in_division - k,
							total_work - current_work - current_work_left,
							num_procs_right, first_proc+num_procs_left,
							first_cell_index+k, proc_index );
				} else {
					right_failure = 0;
					proc_index[first_proc+num_procs_in_division] =
						first_cell_index+num_root_cells_in_division;
				}

				if ( left_failure && right_failure ) {
					return divide_list_linear( global_work, 
							constrained_quantities,
							per_proc_constraints, 
							num_root_cells_in_division,
							total_work, num_procs_in_division, 
							first_proc, first_cell_index,
							proc_index );
				} else if ( left_failure || right_failure ) {
					current_work += global_work[k];
				} else {
					return 0;
				}
			}

			return -1;
		} else {
			return 0;
		}
	}

	cart_debug("i = %d, j = %d, num_procs = %d, %d, %d, num_root_cells = %d",
		i, j, num_procs_in_division, num_procs_left, num_procs_right, 
		num_root_cells_in_division );

	return -1;
}

int divide_list_linear( float *global_work, int *constrained_quantities,
                int *per_proc_constraints, 
		int num_root_cells_in_division, double total_work,
                int num_procs_in_division, int first_proc,
                int first_cell_index, int *proc_index ) {

	int i, c, c2;
	int proc;
	int index;
	double local_work;
	double ideal_work_per_proc;
	int num_procs_remaining;
	long total_constraints[num_constraints];
	long sum_constraints[num_constraints];

	cart_debug("dividing linearly for %u procs, first_proc = %u for %u cells starting at %d", 
			num_procs_in_division, first_proc, num_root_cells_in_division, first_cell_index );

	ideal_work_per_proc = total_work / (double)num_procs_in_division;

	index = 0;
	num_procs_remaining = num_procs_in_division-1;
	proc_index[first_proc] = first_cell_index;

	for ( c = 0; c < num_constraints; c++ ) {
		total_constraints[c] = 0;
	}

	for ( i = 0; i < num_root_cells_in_division; i++ ) {
		for ( c = 0; c < num_constraints; c++ ) {
			total_constraints[c] += constrained_quantities[num_constraints*i+c];
		}
	}

	for ( proc = first_proc+1; proc < num_procs_in_division + first_proc; proc++ ) {
		/* assign a minimum number of cells so that the remaining processors
		 * can satisfy the constraints (cells, particles, etc) */

		/* start by assigning 1 root cell */
		for ( c = 0; c < num_constraints; c++ ) {
			sum_constraints[c] = constrained_quantities[num_constraints*index+c];
			total_constraints[c2] -= constrained_quantities[num_constraints*index+c];
		}
		index++;
	
		for ( c = 0; c < num_constraints; c++ ) {
			while ( total_constraints[c] >= num_procs_remaining*per_proc_constraints[c]  ) {
				for ( c2 = 0; c2 < num_constraints; c2++ ) {
					total_constraints[c2] -= constrained_quantities[num_constraints*index+c2];
					sum_constraints[c2] += constrained_quantities[num_constraints*index+c2];

					if ( sum_constraints[c2] > per_proc_constraints[c2] ) {
						/* problem allocating! */
						cart_debug("failed allocating space linearly!");
						return -1;
					}
				}
				index++;
			}
		}
		
		/* ensure maximum of reduced_root_cells, leaving at least one root cell
		 * for each processor left to come, and not assigning more cells than processor
		 * can hold (max_cells_allowed), as well as allowing to go above ideal_work_per_proc
		 * if rest of procs can't hold the cells remaining
		 */
		while (         index < num_root_cells_in_division - num_procs_remaining &&
				local_work + 0.5*global_work[index] < ideal_work_per_proc ) {

			for ( c = 0; c < num_constraints; c++ ) {
				sum_constraints[c] += constrained_quantities[num_constraints*index+c];

				if ( sum_constraints[c] >= per_proc_constraints[c] ) {
					break;
				}
			}

			if ( c < num_constraints ) {
				break;
			} else {
				for ( c = 0; c < num_constraints; c++ ) {
					total_constraints[c] -= constrained_quantities[num_constraints*index+c];
				}

				local_work += global_work[index];
				index++;
			}
		}

		proc_index[proc] = index+first_cell_index;

		/* continually re-evaluate ideal work given what work is left */
		num_procs_remaining--;
		total_work -= local_work;

		if ( num_procs_remaining > 0 ) {
			ideal_work_per_proc = total_work / (double)(num_procs_remaining);
		}
	}

	proc_index[first_proc+num_procs_in_division] = first_cell_index + num_root_cells_in_division;

	return 0;
}

void load_balance_entire_volume( float *global_work, 
		int *constrained_quantities,
		int *new_proc_sfc_index ) {

	int i, j, c;
	int ret;
	double total_work, current_work;
	double avg_work;
	long sum_constraints[num_constraints];
	int per_proc_constraints[num_constraints];

	char filename[256];
	int foo;
	FILE *output;

#ifdef SAVE_LOAD_BALANCE_PARTITION
	FILE *partition;
#endif

	per_proc_constraints[0] = (1.0-est_buffer_fraction)*(double)num_cells;

#ifdef PARTICLES
	per_proc_constraints[1] = 0.9*num_particles;
#endif /* PARTICLES */

	/* compute total work */
	total_work = 0.0;
	for ( i = 0; i < num_root_cells; i++ ) {
		total_work += global_work[i];
	}

	avg_work = 0.5 * total_work / (double)num_root_cells;
	cart_debug("avg_work = %e", avg_work );

	for ( c = 0; c < num_constraints; c++ ) {
		cart_debug("per_proc_constraints[%u] = %d", c, per_proc_constraints[c] );
	}

	/* find large blocks of low-work (unrefined) cells */
/*
	i = 0;
	count = 0;
	last = -1;
	num_blocks = 0;
	num_reserved_procs = 0;
	while ( i < num_root_cells ) {
		if ( global_work[i] <= avg_work ) {
			count++;
			for ( c = 0; c < num_constraints; c++ ) {
				sum_constraints[c] += constrained_quantities[num_constraints*i+c];
			}
		} else if ( count > 20 ) {
			j = i;
			while ( j < num_root_cells && global_work[j] > avg_work ) {
				j++;
			}

			proc_fraction = 0.0;

			for ( c = 0; c < num_constraints; c++ ) {
				proc_fraction = max( proc_fraction, 
					(float)sum_constraints[c] / (float)per_proc_constraints[c] );
				sum_constraints[c] = 0;
			}

			if ( proc_fraction > 0.2 ) {
				for ( j = last+1; j < i; j++ ) {
					global_work[j] = -global_work[j];
				}
			
				num_blocks++;
				num_reserved_procs += (int)ceil(proc_fraction);
			}

			count = 0;
			last = i;
		} else {
			last = i;
		}

		i++;
	}

	num_reserved_procs = min( num_reserved_procs, num_procs-1 );
	cart_debug("num_reserved_procs = %d", num_reserved_procs );
	
	total_work = 0.0;
	for ( i = 0; i < num_root_cells; i++ ) {
		if ( global_work[i] > 0.0 ) {
			total_work += global_work[i];
		}
	}

	new_proc_sfc_index[0] = 0;
	ideal_work_per_proc = total_work / (double)(num_procs - num_reserved_procs);

	cart_debug("total_work = %e, num_reserved_procs = %u, ideal_work_per_proc = %e", total_work,
		num_reserved_procs, ideal_work_per_proc );

	if ( num_blocks > 0 && num_reserved_procs < num_procs ) {
		last = 0;
		i = 0;
		flag = ( global_work[0] > 0 ) ? 1 : 0;
		current_proc = 0;
		num_procs_in_block = 0;
		while ( i < num_root_cells ) {
			current_work = 0.0;

			if ( flag ) {
				while ( i < num_root_cells && global_work[i] > 0 ) {
					current_work += global_work[i];
					i++;
				}

				num_root_cells_in_block = i - last;
				num_procs_in_block = max( 1, (int)floor(current_work / ideal_work_per_proc) );

				num_reserved_procs += num_procs_in_block;
				total_work -= current_work;
				ideal_work_per_proc = total_work / (double)(num_procs - num_reserved_procs);
			} else {
				while ( i < num_root_cells && global_work[i] < 0.0 ) {
					global_work[i] = -global_work[i];
					current_work += global_work[i];

					for ( c = 0; c < num_constraints; c++ ) {
						sum_constraints[c] += constrained_quantities[num_constraints*i+c];	
					}

					i++;
				}

				num_root_cells_in_block = i - last;

				num_procs_in_block = 0;
				for ( c = 0; c < num_constraints; c++ ) {
					num_procs_in_block = max( num_procs_in_block,
						(int)ceil((float)sum_constraints[c]/(float)per_proc_constraints[c]) );
				}
			}

			num_procs_in_block = min( num_procs_in_block, num_procs - current_proc );

			if ( i < num_root_cells ) {
				flag = ( global_work[i] < 0 ) ? 0 : 1;

				if ( current_proc + num_procs_in_block == num_procs ) {
					for ( ; i < num_root_cells; i++ ) {
						current_work += global_work[i];
					}

					num_root_cells_in_block = num_root_cells - last;
				}
			} else {
				num_procs_in_block = num_procs - current_proc;
			}

			new_proc_sfc_index[current_proc] = last;

			do {
				cart_debug("dividing block starting with current_proc = %d", current_proc );
				ret = divide_list_recursive( &global_work[last], 
						&constrained_quantities[num_constraints*last],
						per_proc_constraints,
						num_root_cells_in_block, current_work, num_procs_in_block,
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
*/
		ret = divide_list_recursive( global_work, 
			constrained_quantities,
			per_proc_constraints,
			num_root_cells, total_work, 
			num_procs, MASTER_NODE, 0, new_proc_sfc_index );
	/* }  */

	if ( ret == -1 ) {
		sprintf( filename, "%s/load_balance.dat", output_directory );
		output = fopen( filename, "w" );

		foo = num_root_cells;
		fwrite( &foo, sizeof(int), 1, output );
		foo = num_constraints;
		fwrite( &foo, sizeof(int), 1, output );
		fwrite( per_proc_constraints, sizeof(int), num_constraints, output );
		fwrite( global_work, sizeof(float), num_root_cells, output );
		fwrite( constrained_quantities, sizeof(int), (long)num_constraints*(long)num_root_cells, output );
		fclose(output);

		cart_error("Unable to find proper load balancing division");

#ifdef SAVE_LOAD_BALANCE_PARTITION
	} else {
		sprintf( filename, "%s/partition.dat", output_directory );
		partition = fopen( filename, "w" );
		cart_assert( partition != NULL );

        fprintf( partition, "%u %u\n", &num_procs, &num_grid );

		for ( i = 0; i < num_procs+1; i++ ) {
			fprintf( partition, "%u\n", &new_proc_sfc_index[i] );
		}

        fclose( partition );
#endif
	}

	total_work = 0.0;
	for ( i = 0; i < num_procs; i++ ) {
		current_work = 0.0;

		for ( c = 0; c < num_constraints; c++ ) {
			sum_constraints[c] = 0;
		}

		for ( j = new_proc_sfc_index[i]; j < new_proc_sfc_index[i+1]; j++ ) {
			current_work += global_work[j];

			for ( c = 0; c < num_constraints; c++ ) {
				sum_constraints[c] += constrained_quantities[num_constraints*j+c];
			}
		}

		total_work += current_work;
		cart_debug("new_proc_sfc_index[%u] = %u, %u, %u, work = %e %d", i, 
				new_proc_sfc_index[i], new_proc_sfc_index[i+1], 
				new_proc_sfc_index[i+1]-new_proc_sfc_index[i],
				current_work, sum_constraints[0] );
	}

	cart_debug("total work = %e", total_work );

	/* sanity tests */
	for ( i = 0; i < num_procs; i++ ) {
		cart_assert( new_proc_sfc_index[i] < new_proc_sfc_index[i+1] );
		cart_assert( new_proc_sfc_index[i+1] <= num_root_cells );
		cart_assert( new_proc_sfc_index[i+1] - new_proc_sfc_index[i] < num_cells );
	}
}

void load_balance() {
	int i;
	int level;
	int root;
	int num_parts;
	float level_cost;
	int proc;
	int first_oct, old_first_oct;
	int new_oct;
	int coords[nDim];
	int receive_counts[MAX_PROCS];
	int receive_displacements[MAX_PROCS];
	int old_proc_sfc_index[MAX_PROCS+1];
	int new_proc_sfc_index[MAX_PROCS+1];

	int ipart;
#ifdef PARTICLES
	int num_parts_to_send[MAX_PROCS];
	int particle_list_to_send[MAX_PROCS];
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	int tracer;
	int tracer_list_to_send[MAX_PROCS];
	int num_tracers_to_send[MAX_PROCS];
#endif /* HYDRO_TRACERS */

	float *local_work;
	int *local_constraints;
	float *global_work;
	int *global_constraints;
	int num_cells_moved;
	int root_cell_shift;
	int icell, ioct;
	int num_level_cells;
	int *level_cells;
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
	start_time( COMMUNICATION_TIMER );

	/* allocate space to store the work estimators */
	local_work = cart_alloc(float, num_cells_per_level[min_level] );
	local_constraints = cart_alloc(int, num_constraints*num_cells_per_level[min_level] );

	#pragma omp parallel for
	for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
		local_work[i] = 0.0;
	}

	#pragma omp parallel for
	for ( i = 0; i < num_constraints*num_cells_per_level[min_level]; i++ ) {
		local_constraints[i] = 0;
	}

	/* compute work and constraints per root cell */
	for ( level = min_level; level <= max_level; level++ ) {
		/* give more weight to higher levels of refinement */
		level_cost = (float)(1<<level); 

		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		#pragma omp parallel for private(icell,root,num_parts,ipart)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			root = cell_parent_root_cell( icell );

			local_work[root] += cost_per_cell * level_cost;
			local_constraints[num_constraints*root]++;

#ifdef PARTICLES
			/* count number of particles in cell */
			num_parts = 0;
			ipart = cell_particle_list[icell];
			while ( ipart != NULL_PARTICLE ) {
				num_parts++;
				ipart = particle_list_next[ipart];				
			}

			local_work[root] += cost_per_particle * (float)num_parts * level_cost;

			local_constraints[num_constraints*root+1] += num_parts;
#endif /* PARTICLES */
		}

		cart_free( level_cells );
	}

	if ( local_proc_id == MASTER_NODE ) {
		global_work = cart_alloc(float, num_root_cells );
		global_constraints = cart_alloc(int, num_constraints * num_root_cells );

		for ( i = 0; i < num_procs; i++ ) {
			receive_counts[i] = proc_sfc_index[i+1] - proc_sfc_index[i];
		}
	}

	/* gather all work to MASTER_NODE */
	start_time( LOAD_BALANCE_COMMUNICATION_TIMER );
	MPI_Gatherv( local_work, num_cells_per_level[min_level], MPI_FLOAT,
		global_work, receive_counts, proc_sfc_index, MPI_FLOAT,
		MASTER_NODE, MPI_COMM_WORLD );
	end_time( LOAD_BALANCE_COMMUNICATION_TIMER );

	if ( local_proc_id == MASTER_NODE ) {
		for ( i = 0; i < num_procs; i++ ) {
			receive_counts[i] *= num_constraints;
			receive_displacements[i] = num_constraints*proc_sfc_index[i];
		}
	}

	start_time( LOAD_BALANCE_COMMUNICATION_TIMER );
	MPI_Gatherv( local_constraints, num_constraints*num_cells_per_level[min_level], MPI_INT,
		global_constraints, receive_counts, receive_displacements, MPI_INT,
		MASTER_NODE, MPI_COMM_WORLD );
	end_time( LOAD_BALANCE_COMMUNICATION_TIMER );

	cart_free( local_work );
	cart_free( local_constraints );

	/* MASTER_NODE determines new splitting of SFC curve */
	if ( local_proc_id == MASTER_NODE ) {
		load_balance_entire_volume( global_work, global_constraints, new_proc_sfc_index );

		cart_free( global_work );
		cart_free( global_constraints );
	}

	/* tell other processors which cells to expect */
	start_time( LOAD_BALANCE_COMMUNICATION_TIMER );	
	MPI_Bcast( &new_proc_sfc_index, num_procs+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	end_time( LOAD_BALANCE_COMMUNICATION_TIMER );

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
	if ( new_proc_sfc_index[local_proc_id] <= proc_sfc_index[local_proc_id] ) {
		if ( new_proc_sfc_index[local_proc_id+1] < proc_sfc_index[local_proc_id+1] ) {
			pack_add_root_trees( transfer_cells, new_proc_sfc_index, 
				max( proc_sfc_index[local_proc_id], new_proc_sfc_index[local_proc_id+1] ),
				proc_sfc_index[local_proc_id+1] );
		}
	} else {
		pack_add_root_trees( transfer_cells, new_proc_sfc_index, 
			proc_sfc_index[local_proc_id], 
			min( new_proc_sfc_index[local_proc_id], proc_sfc_index[local_proc_id+1] ) );

		if ( new_proc_sfc_index[local_proc_id+1] < proc_sfc_index[local_proc_id+1] ) {
			pack_add_root_trees( transfer_cells, new_proc_sfc_index, 
				new_proc_sfc_index[local_proc_id+1],
				proc_sfc_index[local_proc_id+1] );
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
	for ( i = 0; i < num_procs+1; i++ ) {
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
				cart_assert( coords[i] >= 0 && coords[i] < num_grid );
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

	end_time( COMMUNICATION_TIMER );
	end_time( LOAD_BALANCE_TIMER );
}
