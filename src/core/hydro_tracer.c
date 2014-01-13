#include "config.h"
#if defined(HYDRO) && defined(HYDRO_TRACERS)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "io.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "timing.h"
#include "sfc.h"
#include "times.h"
#include "tree.h"
#include "units.h"


double tracer_x[num_tracers][nDim];
tracerid_t tracer_id[num_tracers];
int tracer_list_next[num_tracers];
int tracer_list_prev[num_tracers];

int cell_tracer_list[num_cells];

int num_local_tracers = 0;
tracerid_t num_tracers_total = 0;

int next_free_tracer = 0;
int free_tracer_list = NULL_TRACER;

int tracer_list_enabled = 0;

#if defined(ENRICHMENT) && defined(ENRICHMENT_SNIa)
int num_hydro_vars_traced = 4;
int hydro_vars_traced[] = {
	HVAR_GAS_DENSITY,
	HVAR_INTERNAL_ENERGY,
	//HVAR_GAS_ENERGY,
	//HVAR_PRESSURE,
    //HVAR_ELECTRON_INTERNAL_ENERGY,
	HVAR_METAL_DENSITY_II,
	HVAR_METAL_DENSITY_Ia
};
char *hydro_vars_traced_labels[] = {
	"HVAR_GAS_DENSITY",
	"HVAR_INTERNAL_ENERGY",    // internal energy
	//"HVAR_GAS_ENERGY",       // total energy
	//"HVAR_PRESSURE",
    //"HVAR_ELECTRON_INTERNAL_ENERGY",
	"HVAR_METAL_DENSITY_II",
	"HVAR_METAL_DENSITY_Ia"
};
#else  /* ENRICHMENT && ENRICHMENT_SNIa */
int num_hydro_vars_traced = 2;
int hydro_vars_traced[] = {
    HVAR_GAS_DENSITY,
    HVAR_INTERNAL_ENERGY
    //HVAR_GAS_ENERGY,
    //HVAR_PRESSURE,
    //HVAR_ELECTRON_INTERNAL_ENERGY
};
char *hydro_vars_traced_labels[] = {
    "HVAR_GAS_DENSITY",
    "HVAR_INTERNAL_ENERGY"    // internal energy
    //"HVAR_GAS_ENERGY",       // total energy
    //"HVAR_PRESSURE",
    //"HVAR_ELECTRON_INTERNAL_ENERGY"
};
#endif /* ENRICHMENT && ENRICHMENT_SNIa */


void init_hydro_tracers() {
	int i;

#pragma omp parallel for default(none), private(i), shared(tracer_id,tracer_list_next,tracer_list_prev)
	for ( i = 0; i < num_tracers; i++ ) {
		tracer_id[i] = NULL_TRACER;
		tracer_list_next[i] = NULL_TRACER;
		tracer_list_prev[i] = NULL_TRACER;
	}

#pragma omp parallel for default(none), private(i), shared(cell_tracer_list,size_cell_array)
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
	double pos[nDim];

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

	MPI_Allgather( &num_leafs, 1, MPI_INT, proc_num_leafs, 1, MPI_INT, mpi.comm.run );

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
				cell_center_position( icell, pos );
				for ( j = 0; j < nDim; j++ ) {
					tracer_x[tracer][j] = pos[j];
				}
				insert_tracer( icell, tracer );
				id++;
			}
		}

		cart_free( level_cells );
	}

	tracer_list_enabled = 1;
}

#ifdef PARTICLES
void set_hydro_tracers_to_particles() {
	int i;
	int ipart;
	int icell;
	int tracer;
	tracerid_t tmp;

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && particle_id[ipart] < particle_species_indices[1] ) {
			if ( particle_id[ipart] > TRACERID_MAX ) {
				cart_error("Integer overflow detecdted in particle to tracer ids. Try compiling with 64-bit tracer ids!");
			}
			tracer = tracer_alloc( particle_id[ipart] );

			for ( i = 0; i < nDim; i++ ) {
				tracer_x[tracer][i] = particle_x[ipart][i];
			}

			icell = cell_find_position( tracer_x[tracer] );
			cart_assert( icell != NULL_OCT );
			insert_tracer( icell, tracer );
		}
	}

	tmp = num_local_tracers;
	MPI_Allreduce( &tmp, &num_tracers_total, 1, MPI_TRACERID_T, MPI_SUM, mpi.comm.run );

	cart_debug("num_local_tracers = %u", num_local_tracers );
	cart_debug("num_tracers_total = %ld", num_tracers_total );
	tracer_list_enabled = 1;
}
#endif /* PARTICLES */

void update_tracer_list( int level ) {
	int i, k;
	int tracer;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	int num_tracers_to_send[MAX_PROCS];
	int tracer_list_to_send[MAX_PROCS];
	int tracer2, new_cell;
	int last;
	int proc;
	int collect_level;
	int sfc;

	start_time( WORK_TIMER );

	/* now move tracers from one cell list to another */
	for ( i = 0; i < num_procs; i++ ) {
		num_tracers_to_send[i] = 0;
		tracer_list_to_send[i] = NULL_TRACER;
	}

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( k = 0; k < num_level_cells; k++ ) {
		iter_cell = level_cells[k];
		tracer = cell_tracer_list[iter_cell];

		while ( tracer != NULL_TRACER ) {
			tracer2 = tracer_list_next[tracer];

			sfc = sfc_index_position( tracer_x[tracer] );
			proc = processor_owner( sfc );

			if ( proc == local_proc_id ) {
				new_cell = cell_find_position_sfc( sfc, tracer_x[tracer] );

				if ( new_cell != iter_cell ) {
					delete_tracer( iter_cell, tracer );
					insert_tracer( new_cell, tracer );
				}
			} else if ( proc == -1 ) {
				cart_error("Unable to find processor owner for tracer %ld", tracer_id[tracer] );
			} else {
				delete_tracer( iter_cell, tracer );
				tracer_list_next[tracer] = tracer_list_to_send[proc];
				tracer_list_to_send[proc] = tracer;
				num_tracers_to_send[proc]++;
			}

			tracer = tracer2;
		}
	}

	cart_free( level_cells );
	end_time( WORK_TIMER );

	start_time( COMMUNICATION_TIMER );
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
	tracerid_t *send_id;
	tracerid_t *recv_id[MAX_PROCS];
	double *send_tracers;
	double *recv_tracers[MAX_PROCS];
	int num_pages_received;

	MPI_Request *send_requests;
	MPI_Request recv_id_requests[MAX_PROCS];
	MPI_Request recv_tracer_requests[MAX_PROCS];
	MPI_Status status;

	page_size = MIN(65536/num_procs, 1024);
	tracers_page_size = nDim*page_size;

	/* set up receives */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( ( trade_level == -1 && proc != local_proc_id ) ||
				( trade_level != -1 &&
				( num_remote_buffers[trade_level][proc] > 0 || num_local_buffers[trade_level][proc] > 0 ) ) ) {
			recv_id[proc] = cart_alloc(tracerid_t, page_size );
			recv_tracers[proc] = cart_alloc(double, tracers_page_size );

			MPI_Irecv( recv_id[proc], page_size, MPI_TRACERID_T, proc, 0,
					mpi.comm.run, &recv_id_requests[proc] );
			MPI_Irecv( recv_tracers[proc], tracers_page_size, MPI_DOUBLE,
					proc, 0, mpi.comm.run, &recv_tracer_requests[proc] );

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
	send_id = cart_alloc(tracerid_t, num_pages_to_send * page_size );
	send_tracers = cart_alloc(double, num_pages_to_send * tracers_page_size );
	send_requests = cart_alloc(MPI_Request, 2*num_pages_to_send );

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
						MPI_Isend( &send_id[num_pages_sent*page_size], page_size, MPI_TRACERID_T, proc,
							proc_pages_sent, mpi.comm.run, &send_requests[num_send_requests++] );
						MPI_Isend( &send_tracers[num_pages_sent*tracers_page_size], tracers_page_size,
							MPI_DOUBLE, proc, proc_pages_sent, mpi.comm.run,
							&send_requests[num_send_requests++] );

						proc_pages_sent++;
						num_pages_sent++;
						id_count = 0;
						tracer_count = 0;
					}
				}

				tracer_list_free( tracer_list_to_send[proc] );

				MPI_Isend( &send_id[num_pages_sent*page_size], id_count, MPI_TRACERID_T, proc, proc_pages_sent,
					mpi.comm.run, &send_requests[num_send_requests++] );
				MPI_Isend( &send_tracers[num_pages_sent*tracers_page_size], tracer_count, MPI_DOUBLE, proc,
					proc_pages_sent, mpi.comm.run, &send_requests[num_send_requests++] );

				num_pages_sent++;
				proc_pages_sent++;

			} else {
				/* send a single empty page */
				MPI_Isend( &send_id[page_size*num_pages_sent], 0, MPI_TRACERID_T, proc, 0, mpi.comm.run,
					&send_requests[num_send_requests++] );
				MPI_Isend( &send_tracers[tracers_page_size*num_pages_sent], 0, MPI_DOUBLE, proc, 0,
					mpi.comm.run, &send_requests[num_send_requests++] );

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
			MPI_Get_count( &status, MPI_TRACERID_T, &id_count );

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
				MPI_Irecv( recv_id[proc], page_size, MPI_TRACERID_T, proc, page_count[proc],
						mpi.comm.run, &recv_id_requests[proc] );
				MPI_Irecv( recv_tracers[proc], tracers_page_size, MPI_DOUBLE, proc,
						page_count[proc], mpi.comm.run, &recv_tracer_requests[proc] );
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

int tracer_alloc( tracerid_t id ) {
	int tracer;

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
	double pos[nDim];

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( cell_is_refined(icell) );

	tracer = cell_tracer_list[icell];
	cell_tracer_list[icell] = NULL_TRACER;

	cell_center_position( icell, pos );

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
	cart_assert( cell_is_local( icell ) );
	cart_assert( cell_is_leaf( icell ) );
	cart_assert( tracer >= 0 && tracer < num_tracers );
	cart_assert( tracer_id[tracer] != NULL_TRACER );

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

#endif /* HYDRO && HYDRO_TRACERS */
