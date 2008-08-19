#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "particle.h"
#include "timestep.h"
#include "iterators.h"
#include "constants.h"
#include "parallel.h"
#include "auxiliary.h"
#include "timing.h"
#include "io.h"
#include "sfc.h"
#include "units.h"

#ifdef HYDRO_TRACERS

#ifndef HYDRO
#error "Hydro tracer particles requires HYDRO to be defined!"
#endif /* HYDRO */

double tracer_x[num_tracers][nDim];
int tracer_id[num_tracers];
int tracer_list_next[num_tracers];
int tracer_list_prev[num_tracers];

int cell_tracer_list[num_cells];

int num_tracer_row    = 256;
int num_local_tracers = 0;
int num_tracers_total = 0;

int next_free_tracer = 0;
int free_tracer_list = NULL_TRACER;

int tracer_list_enabled = 0;

#ifdef ENRICH
int num_hydro_vars_traced = 4;
int hydro_vars_traced[] = {		HVAR_GAS_DENSITY,
					HVAR_INTERNAL_ENERGY,
					HVAR_METALLICITY_II,
					HVAR_METALLICITY_Ia };
char *hydro_vars_traced_labels[] = { 	"density",
					"internal energy",
					"SNII metallicity",
					"SNIa metallicity" };
#else 
int num_hydro_vars_traced = 3;
int hydro_vars_traced[] = {      	HVAR_GAS_DENSITY,
					HVAR_GAS_ENERGY,
                                        HVAR_INTERNAL_ENERGY };
char *hydro_vars_traced_labels[] = {	"density",
					"total energy",
					"internal energy" };
#endif /* ENRICH */


void init_hydro_tracers() { 
	int i;

	for ( i = 0; i < num_tracers; i++ ) {
		tracer_id[i] = NULL_TRACER;
	}

#pragma omp parallel for default(none), private(i), shared(cell_tracer_list)
	for ( i = 0; i < num_cells; i++ ) {
		cell_tracer_list[i] = NULL_TRACER;
	}

	num_local_tracers = 0;
	num_tracers_total = 0;
	next_free_tracer = 0;
	free_tracer_list = NULL_TRACER;

	tracer_list_enabled = 0;
}

void set_hydro_tracers( int min_tracer_level ) {
	int i, j;
	int level;
	int num_level_cells;
	int *level_cells;
	int icell;
	int tracer;
	unsigned int id;
	int num_leafs;
	int proc;
	int proc_num_leafs[MAX_PROCS];
	float pos[nDim];

	num_leafs = 0;
	for ( level = min_tracer_level; level < max_level; level++ ) {
		num_leafs += num_cells_per_level[level] - num_cells_per_level[level+1] / num_children;
	}

	if ( max_level >= min_tracer_level ) {
		num_leafs += num_cells_per_level[max_level];
	}

	cart_debug("num_leafs = %u", num_leafs );

	if ( num_leafs > num_tracers ) {
		cart_error("num_tracers < num_leafs in set_hydro_tracers, increase num_tracers" );
	}

	MPI_Allgather( &num_leafs, 1, MPI_INT, proc_num_leafs, 1, MPI_INT, MPI_COMM_WORLD );

	id = 0;
	for ( proc = 0; proc < local_proc_id; proc++ ) {
		id += proc_num_leafs[proc];
	}

	num_tracers_total = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		num_tracers_total += proc_num_leafs[proc];
	}
	
	/* put a tracer particle at the center of every leaf cell */
	for ( level = min_tracer_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf(icell) ) {
				tracer = tracer_alloc( id );
				cell_position( icell, pos );
				for ( j = 0; j < nDim; j++ ) {
					tracer_x[tracer][j] = pos[j];
				}
				insert_tracer( icell, tracer );
				id++;
			}
		}
	}

	tracer_list_enabled = 1;
}

#ifdef PARTICLES
void set_hydro_tracers_to_particles() { 
	int i;
	int ipart;
	int icell;
	int tracer;

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && particle_id[ipart] < particle_species_indices[1] ) {
			tracer = tracer_alloc( particle_id[ipart] );

			for ( i = 0; i < nDim; i++ ) {
				tracer_x[tracer][i] = particle_x[ipart][i];
			}

			icell = cell_find_position( tracer_x[tracer] );
			insert_tracer( icell, tracer );
		}
	}

	MPI_Allreduce( &num_local_tracers, &num_tracers_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

	cart_debug("num_local_tracers = %u", num_local_tracers );
	cart_debug("num_tracers_total = %u", num_tracers_total );
	tracer_list_enabled = 1;
}
#endif /* PARTICLES */

void move_hydro_tracers( int level ) {
	int i, j;
	int tracer;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	double vdt[nDim];
	int icell, icell_orig;
	int level1;
	int child;
	float pos[nDim];
	int found;
	int c[num_children];
	double diff1, diff2, diff3;
	double pt3, pd3;
	double t1,t2,t3,d1,d2,d3;
	double t2t1, t2d1, d2t1, d2d1;
	double t3t2t1, t3t2d1, t3d2t1, t3d2d1;
	double d3t2t1, d3t2d1, d3d2t1, d3d2d1;

	start_time( WORK_TIMER ); 

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];

#ifdef HYDRO_TRACERS_NGP
		for ( j = 0; j < nDim; j++ ) {
			vdt[j] = cell_momentum(iter_cell,j)/cell_gas_density(iter_cell) * dtl[level];
		}
#endif /* HYDRO_TRACERS_NGP */

		tracer = cell_tracer_list[iter_cell];
		while ( tracer != NULL_TRACER ) {
			cart_assert( tracer >= 0 && tracer < num_tracers );

#ifndef HYDRO_TRACERS_NGP
			icell = iter_cell;
			level1 = level;

			do {
				found = 1;
				icell_orig = icell;
				cart_assert( icell != NULL_OCT );

				cell_position( icell, pos );

				/* find lower leftmost cell */
				child = 0;
				for ( j = 0; j < nDim; j++ ) {
					if ( tracer_x[tracer][j] >= pos[j] ) {
						child += (1<<j);
					}
				}

				cart_assert( child >= 0 && child < num_children );

				for ( j = 0; j < nDim; j++ ) {
					if ( neighbor_moves[child][j] == -1 ) {
						break;
					} else {
						icell = cell_neighbor(icell, neighbor_moves[child][j] );
						cart_assert( icell != NULL_OCT );

						if ( cell_level(icell) != level1 ) {
							icell = cell_parent_cell(icell_orig);
							cart_assert( icell != NULL_OCT );
							level1 = level1 - 1;
							found = 0;
							break;
						}
					}
				}

				if ( found ) {
					c[0] = icell;
					c[1] = cell_neighbor(icell,1);
					c[2] = cell_neighbor(icell,3);
					c[3] = cell_neighbor(c[1],3);
					c[4] = cell_neighbor(icell,5);
					c[5] = cell_neighbor(c[1],5);
					c[6] = cell_neighbor(c[2],5);
					c[7] = cell_neighbor(c[3],5);

					for ( j = 1; j < num_children; j++ ) {
						if ( cell_level(c[j]) != level1 ) {
							icell = cell_parent_cell(icell_orig);
							level1 = level1 - 1;
							cart_assert( icell != NULL_OCT );
							found = 0;
							break;
						}
					}
				}
			} while ( !found );

			cell_position( c[0], pos );

			/* now we have the level on which this particle will move */
			diff1 = pos[0] - tracer_x[tracer][0];
			if ( fabs(diff1) > (double)(num_grid/2) ) {
				if ( diff1 > 0.0 ) {
					diff1 -= (double)(num_grid);
				} else {
					diff1 += (double)(num_grid);
				}
			}
			d1 = fabs(diff1) * cell_size_inverse[level1];
			cart_assert( d1 >= 0.0 && d1 <= 1.0 );

			diff2 = pos[1] - tracer_x[tracer][1];
			if ( fabs(diff2) > (double)(num_grid/2) ) {
				if ( diff2 > 0.0 ) {
					diff2 -= (double)(num_grid);
				} else {
					diff2 += (double)(num_grid);
				}
			}
			d2 = fabs(diff2) * cell_size_inverse[level1];

			diff3 = pos[2] - tracer_x[tracer][2];
			if ( fabs(diff3) > (double)(num_grid/2) ) {
				if ( diff3 > 0.0 ) {
					diff3 -= (double)(num_grid);
				} else {
					diff3 += (double)(num_grid);
				}
			}
			d3 = fabs(diff3) * cell_size_inverse[level1];

			cart_assert( d1 >= 0.0 && d1 <= 1.0 );
			cart_assert( d2 >= 0.0 && d2 <= 1.0 );
			cart_assert( d3 >= 0.0 && d3 <= 1.0 );

			t1   = 1.0 - d1;
			t2   = 1.0 - d2;
			t3   = 1.0 - d3;

			cart_assert( t1 >= 0.0 && t1 <= 1.0 );
			cart_assert( t2 >= 0.0 && t2 <= 1.0 );
			cart_assert( t3 >= 0.0 && t3 <= 1.0 );

			t2t1 = t2 * t1;
			t2d1 = t2 * d1;
			d2t1 = d2 * t1;
			d2d1 = d2 * d1;

			pt3 = t3*dtl[level];
			pd3 = d3*dtl[level];

			t3t2t1 = pt3 * t2t1;
			t3t2d1 = pt3 * t2d1;
			t3d2t1 = pt3 * d2t1;
			t3d2d1 = pt3 * d2d1;
			d3t2t1 = pd3 * t2t1;
			d3t2d1 = pd3 * t2d1;
			d3d2t1 = pd3 * d2t1;
			d3d2d1 = pd3 * d2d1;

			for ( j = 0; j < nDim; j++ ) {
				vdt[j] =t3t2t1 * cell_momentum(c[0], j) / cell_gas_density(c[0]) +
					t3t2d1 * cell_momentum(c[1], j) / cell_gas_density(c[1]) +
					t3d2t1 * cell_momentum(c[2], j) / cell_gas_density(c[2]) +
					t3d2d1 * cell_momentum(c[3], j) / cell_gas_density(c[3]) +
					d3t2t1 * cell_momentum(c[4], j) / cell_gas_density(c[4]) +
					d3t2d1 * cell_momentum(c[5], j) / cell_gas_density(c[5]) +
					d3d2t1 * cell_momentum(c[6], j) / cell_gas_density(c[6]) +
					d3d2d1 * cell_momentum(c[7], j) / cell_gas_density(c[7]);
			}
#endif /* HYDRO_TRACERS_NGP */

			tracer_x[tracer][0] += vdt[0];
			tracer_x[tracer][1] += vdt[1];
			tracer_x[tracer][2] += vdt[2];

			/* enforce periodic boundaries */
			if ( tracer_x[tracer][0] < 0.0 ) {
				tracer_x[tracer][0] += (double)(num_grid);		
			}
				
			if ( tracer_x[tracer][0] >= (double)(num_grid) ) {
				tracer_x[tracer][0] -= (double)(num_grid);
			}

			if ( tracer_x[tracer][1] < 0.0 ) {
				tracer_x[tracer][1] += (double)(num_grid);
			}

			if ( tracer_x[tracer][1] >= (double)(num_grid) ) {
				tracer_x[tracer][1] -= (double)(num_grid);
			}

			if ( tracer_x[tracer][2] < 0.0 ) {
				tracer_x[tracer][2] += (double)(num_grid);
			}

			if ( tracer_x[tracer][2] >= (double)(num_grid) ) {
				tracer_x[tracer][2] -= (double)(num_grid);
			}

			tracer = tracer_list_next[tracer];
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
}

void update_tracer_list( int level ) {
	int i, j, k;
	int tracer, icell;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	int coords[nDim];
	int num_tracers_to_send[MAX_PROCS];
	int tracer_list_to_send[MAX_PROCS];
	int tracer2, new_cell;
	int last;
	int proc;
	int collect_level;

	int min_level_modified, max_level_modified;

	start_time( WORK_TIMER );

	min_level_modified = max_level_modified = level;

        /* now move tracers from one cell list to another */
        for ( i = 0; i < num_procs; i++ ) {
		num_tracers_to_send[i] = 0;
		tracer_list_to_send[i] = NULL_TRACER;
        }

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( k = 0; k < num_level_cells; k++ ) {
		iter_cell = level_cells[k];

		if ( cell_is_leaf(iter_cell) ) {
			tracer = cell_tracer_list[iter_cell];

			while ( tracer != NULL_TRACER ) {
        	                tracer2 = tracer_list_next[tracer];
                	        new_cell = cell_find_position( tracer_x[tracer] );

				if ( new_cell == NULL_OCT ) {
					for ( i = 0; i < nDim; i++ ) {
						coords[i] = (int)(tracer_x[tracer][i]);
					}
					proc = processor_owner( sfc_index( coords ) );

					delete_tracer( iter_cell, tracer );
					tracer_list_next[tracer] = tracer_list_to_send[proc];
					tracer_list_to_send[proc] = tracer;
				} else { 
					if ( new_cell != iter_cell ) {
						if ( cell_level(new_cell) < min_level_modified ) {
							min_level_modified = cell_level( new_cell );
						}
	
						if ( cell_level(new_cell) > max_level_modified ) {
							max_level_modified = cell_level( new_cell );
						}
	
						delete_tracer( iter_cell, tracer );
						insert_tracer( new_cell, tracer );
					}
				}

				tracer = tracer2;
        	        }
	        }
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );

	start_time( COMMUNICATION_TIMER );

	/* now collect tracers which ended up in buffer cells */
	for ( collect_level = min_level_modified; collect_level <= max_level_modified; collect_level++ ) {
		select_level( collect_level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			iter_cell = level_cells[i];

			if ( cell_tracer_list[iter_cell] != NULL_TRACER ) {
				proc = processor_owner( cell_parent_root_sfc( iter_cell ) );

				last = cell_tracer_list[iter_cell];
				num_tracers_to_send[proc]++;

				while ( tracer_list_next[last] != NULL_TRACER ) {
					last = tracer_list_next[last];
					num_tracers_to_send[proc]++;
				}

		                tracer_list_next[last] = tracer_list_to_send[proc];
		                tracer_list_to_send[proc] = cell_tracer_list[iter_cell];
                                                                                                                                                            
                		cell_tracer_list[iter_cell] = NULL_TRACER;
		        }
		}
		cart_free( level_cells );
	}

	trade_tracer_lists( num_tracers_to_send, tracer_list_to_send, level );

	end_time( COMMUNICATION_TIMER );
}

void trade_tracer_lists( int *num_tracers_to_send, int *tracer_list_to_send, int trade_level ) {
	int i, j;
	int id_count, tracer_count;
	int proc, icell, tracer;
	int page_size, tracers_page_size;
	int proc_pages_sent;
	int num_pages_sent;
	int num_pages_to_send;
	int num_send_requests;
	int page_count[MAX_PROCS];
	unsigned int *send_id;
	unsigned int *recv_id[MAX_PROCS];
	double *send_tracers;
	double *recv_tracers[MAX_PROCS];
	int num_pages_received;

	MPI_Request *send_requests;
	MPI_Request recv_id_requests[MAX_PROCS];
	MPI_Request recv_tracer_requests[MAX_PROCS];
	MPI_Status status;

	page_size = (num_tracer_row*num_tracer_row)/num_procs;
	tracers_page_size = nDim*page_size;

	/* set up receives */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( ( trade_level == -1 && proc != local_proc_id ) ||
				( trade_level != -1 &&
				( num_remote_buffers[trade_level][proc] > 0 || num_local_buffers[trade_level][proc] > 0 ) ) ) {
			recv_id[proc] = cart_alloc( page_size * sizeof(unsigned int) );
			recv_tracers[proc] = cart_alloc( tracers_page_size * sizeof(double) );

			MPI_Irecv( recv_id[proc], page_size, MPI_INT, proc, 0, 
					MPI_COMM_WORLD, &recv_id_requests[proc] );
			MPI_Irecv( recv_tracers[proc], tracers_page_size, MPI_DOUBLE, 
					proc, 0, MPI_COMM_WORLD, &recv_tracer_requests[proc] );

			page_count[proc] = 0;
		} else {
			recv_id_requests[proc] = MPI_REQUEST_NULL;
			recv_tracer_requests[proc] = MPI_REQUEST_NULL;
		}
	}

	/* compute number of pages we'll need to send */
	num_pages_to_send = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( ( trade_level == -1 && proc != local_proc_id ) || 
				( trade_level != -1 && 				
				( num_remote_buffers[trade_level][proc] > 0 
				|| num_local_buffers[trade_level][proc] > 0 ) ) ) {
			/* must send at least 1 page */
			if ( num_tracers_to_send[proc] > 0 ) {
				num_pages_to_send += ( num_tracers_to_send[proc] ) / page_size + 1;
			} else {
				num_pages_to_send++;
			}
		}
	}

	/* allocate space for the pages */
	send_id = cart_alloc( num_pages_to_send * page_size * sizeof(int) );
	send_tracers = cart_alloc( num_pages_to_send * tracers_page_size * sizeof(double) );
	send_requests = cart_alloc( 2*num_pages_to_send * sizeof(MPI_Request) );
	
	num_pages_sent = 0;
	num_send_requests = 0;
        for ( proc = 0; proc < num_procs; proc++ ) {
		if ( ( trade_level == -1 && proc != local_proc_id ) || 
				( trade_level != -1 &&  
				( num_remote_buffers[trade_level][proc] > 0 || 
				num_local_buffers[trade_level][proc] > 0 ) ) ) {
			if ( num_tracers_to_send[proc] > 0 ) {
				tracer_count = 0;
				id_count = 0;
				proc_pages_sent = 0;

				tracer = tracer_list_to_send[proc];
	
				while ( tracer != NULL_TRACER ) {
					send_id[page_size*num_pages_sent+id_count++] = tracer_id[tracer];
	
					for ( j = 0; j < nDim; j++ ) {
						send_tracers[tracers_page_size*num_pages_sent+tracer_count++] = 
							tracer_x[tracer][j];
					}
	
					tracer = tracer_list_next[tracer];

					if ( id_count == page_size ) {
						MPI_Isend( &send_id[num_pages_sent*page_size], page_size, MPI_INT, proc, 
							proc_pages_sent, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
						MPI_Isend( &send_tracers[num_pages_sent*tracers_page_size], tracers_page_size, 
							MPI_DOUBLE, proc, proc_pages_sent, MPI_COMM_WORLD, 
							&send_requests[num_send_requests++] );

						proc_pages_sent++;
						num_pages_sent++;
						id_count = 0;
						tracer_count = 0;
					}
				}
	
				tracer_list_free( tracer_list_to_send[proc] );
	
				MPI_Isend( &send_id[num_pages_sent*page_size], id_count, MPI_INT, proc, proc_pages_sent, 
					MPI_COMM_WORLD, &send_requests[num_send_requests++] );
				MPI_Isend( &send_tracers[num_pages_sent*tracers_page_size], tracer_count, MPI_DOUBLE, proc, 
					proc_pages_sent, MPI_COMM_WORLD, &send_requests[num_send_requests++] );

				num_pages_sent++;
				proc_pages_sent++;

			} else {
				/* send a single empty page */
				MPI_Isend( &send_id[page_size*num_pages_sent], 0, MPI_INT, proc, 0, MPI_COMM_WORLD, 
					&send_requests[num_send_requests++] );
				MPI_Isend( &send_tracers[tracers_page_size*num_pages_sent], 0, MPI_DOUBLE, proc, 0, 
					MPI_COMM_WORLD, &send_requests[num_send_requests++] );

				num_pages_sent++;
			}
		} else {
			cart_assert( num_tracers_to_send[proc] == 0 );
		}
	}

	cart_assert( num_pages_sent == num_pages_to_send );

	/* wait for receives and process tracers */
	num_pages_received = 0;
	do {
		MPI_Waitany( num_procs, recv_id_requests, &proc, &status );

		if ( proc != MPI_UNDEFINED ) {
			num_pages_received++;
			MPI_Get_count( &status, MPI_INT, &id_count );

			MPI_Wait( &recv_tracer_requests[proc], MPI_STATUS_IGNORE );

			/* process received page */
			tracer_count = 0;
			for ( i = 0; i < id_count; i++ ) {
				cart_assert( recv_id[proc][i] < num_tracers_total );
				tracer = tracer_alloc( recv_id[proc][i] );

				for ( j = 0; j < nDim; j++ ) {
					tracer_x[tracer][j] = recv_tracers[proc][tracer_count++];
				}

				icell = cell_find_position( tracer_x[tracer] );
				cart_assert( icell != -1 && cell_is_local(icell) );
				insert_tracer( icell, tracer );
			}

			/* if we received a full page, set up to receive a new one */
                        if ( id_count == page_size ) {
                                page_count[proc]++;
				MPI_Irecv( recv_id[proc], page_size, MPI_INT, proc, page_count[proc], 
					MPI_COMM_WORLD, &recv_id_requests[proc] );
				MPI_Irecv( recv_tracers[proc], tracers_page_size, MPI_DOUBLE, proc, 
					page_count[proc], MPI_COMM_WORLD, &recv_tracer_requests[proc] );
                        } else {
				cart_free( recv_id[proc] );
				cart_free( recv_tracers[proc] );
			}
		}
	} while ( proc != MPI_UNDEFINED );

	/* wait for sends to complete */
	MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );

	/* de-allocate send buffers */
	cart_free( send_id );
	cart_free( send_tracers );
	cart_free( send_requests );
}

void build_tracer_list() {
	int i;
	int icell;
	float pos[nDim];

	cart_debug("build_tracer_list()");

	for ( i = 0; i < num_tracers; i++ ) {
		if ( tracer_id[i] != NULL_TRACER ) {
			/* find cell tracer belongs to */
			icell = cell_find_position( tracer_x[i] );
			insert_tracer( icell, i );
		}
	}

	tracer_list_enabled = 1;
}

int tracer_alloc( int id ) {
	int tracer;
	int i;

	if ( free_tracer_list == NULL_TRACER ) {
		if ( num_local_tracers >= num_tracers ) {
			/* generate an error, ran out of tracers */
			cart_error("Ran out of hydro tracer particles, increase num_tracers!");
		}

		tracer = next_free_tracer;
		next_free_tracer++;
	} else {
		tracer = free_tracer_list;
		free_tracer_list = tracer_list_next[free_tracer_list];
	}

	tracer_id[tracer] = id;
	num_local_tracers++;

	return tracer;
}

void tracer_free( int tracer ) {
	cart_assert( tracer >= 0 && tracer < num_tracers );
	cart_assert( tracer_id[tracer] >= 0 );

	tracer_id[tracer] = NULL_TRACER;

	num_local_tracers--;

	tracer_list_next[tracer] = free_tracer_list;
	free_tracer_list = tracer;
}

void tracer_list_free( int ihead ) {
	int last, next;

	cart_assert( ihead != NULL_TRACER );
	last = ihead;

	while ( last != NULL_TRACER ) {
		next = tracer_list_next[last];

		tracer_list_next[last] = free_tracer_list;
		free_tracer_list = last;
		tracer_id[last] = NULL_TRACER;
		num_local_tracers--;

		last = next;
	}
}

void split_tracer_list( int icell ) {
	int i;
	int tracer;
	int next;
	int child;
	float pos[nDim];

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( cell_is_refined(icell) );

	tracer = cell_tracer_list[icell];
	cell_tracer_list[icell] = NULL_TRACER;

	cell_position( icell, pos );

	while ( tracer != NULL_TRACER ) {
		next = tracer_list_next[tracer];

		/* find which child contains this tracer */
		child = 0;
		for ( i = 0; i < nDim; i++ ) {
			if ( tracer_x[tracer][i] >= pos[i] ) {
				child += (1<<i);
			}
		}

		insert_tracer( cell_child( icell, child ), tracer );
		tracer = next;
	}
}

void join_tracer_list( int icell ) {
	int i;
	int tracer;
	int head;

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( cell_is_refined(icell) );
	
	/* find children with tracers */
	tracer = NULL_TRACER;
	for ( i = 0; i < num_children; i++ ) {
		head = cell_tracer_list[ cell_child( icell, i ) ];
		cell_tracer_list[ cell_child( icell, i ) ] = NULL_TRACER;

		if ( head != NULL_TRACER ) {
			if ( tracer == NULL_TRACER ) {
				cell_tracer_list[icell] = head;
			} else {
				tracer_list_next[tracer] = head;
				tracer_list_prev[head] = tracer;
			}

			tracer = head;
			while ( tracer_list_next[tracer] != NULL_TRACER ) {
				tracer = tracer_list_next[tracer];
			}
		}
	}
}

void insert_tracer( int icell, int tracer ) {
	int head;

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( cell_is_leaf( icell ) );
	cart_assert( tracer >= 0 && tracer < num_tracers );

	head = cell_tracer_list[icell];

	tracer_list_prev[tracer] = NULL_TRACER;
	tracer_list_next[tracer] = head;

	if ( head != NULL_TRACER ) {
		cart_assert( tracer_list_prev[head] == NULL_TRACER );
		tracer_list_prev[head] = tracer;
	}
	
	cell_tracer_list[icell] = tracer;
}

void delete_tracer( int icell, int tracer ) {
	int next, prev;

	cart_assert( tracer >= 0 && tracer < num_tracers );

	prev = tracer_list_prev[tracer];
	next = tracer_list_next[tracer];

	if ( prev == NULL_TRACER ) {
		cell_tracer_list[icell] = next;
	} else {
		tracer_list_next[prev] = next;
	}

	if ( next != NULL_TRACER ) {
		tracer_list_prev[next] = prev;
	}
}

#endif /* HYDRO_TRACERS */
