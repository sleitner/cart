#include "config.h"
#ifdef PARTICLES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "density.h"
#include "hydro.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "sfc.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


double particle_t[num_particles];
double particle_dt[num_particles];
double particle_x[num_particles][nDim];
double particle_v[num_particles][nDim];

#ifdef GRAVITY
float particle_pot[num_particles];
#endif /* GRAVITY */

int particle_level[num_particles];
float particle_mass[num_particles];
particleid_t particle_id[num_particles];
int particle_list_next[num_particles];
int particle_list_prev[num_particles];

/* particle species */
int num_particle_species = 0;
float particle_species_mass[MAX_PARTICLE_SPECIES];
particleid_t particle_species_num[MAX_PARTICLE_SPECIES];
particleid_t particle_species_indices[MAX_PARTICLE_SPECIES+1];

/* variables for logging energy */
double tintg = 0.0;
double ekin = 0.0;
double ekin1 = 0.0;
double ekin2 = 0.0;
double au0 = 0.0;
double aeu0 = 0.0;
double ap0 = 0.0;
double ap1 = 0.0;

int cell_particle_list[num_cells];

int num_local_particles = 0;
particleid_t num_particles_total = 0;

int next_free_particle = 0;
int free_particle_list = NULL_PARTICLE;

#ifdef STAR_FORMATION
int next_free_star_particle = 0;
int free_star_particle_list = NULL_PARTICLE;
#endif /* STAR_FORMATION */

int particle_list_enabled = 0;

void init_particles() { 
	int i;

#pragma omp parallel for default(none), private(i), shared(cell_particle_list,size_cell_array)
	for ( i = 0; i < num_cells; i++ ) {
		cell_particle_list[i] = NULL_PARTICLE;
	}

#pragma omp parallel for default(none), private(i), shared(particle_id,particle_level,particle_list_next,particle_list_prev)
	for ( i = 0; i < num_particles; i++ ) {
		particle_id[i] = NULL_PARTICLE;
		particle_level[i] = FREE_PARTICLE_LEVEL;
		particle_list_next[i] = NULL_PARTICLE;
		particle_list_prev[i] = NULL_PARTICLE;
	}

	num_local_particles = 0;
	free_particle_list = NULL_PARTICLE;

#ifdef STAR_FORMATION
	free_star_particle_list = NULL_PARTICLE;
	next_free_star_particle = 0;
	next_free_particle = num_star_particles;
	num_local_star_particles = 0;
#else 
	next_free_particle = 0;
#endif /* STAR_FORMATION */

	particle_list_enabled = 0;
}

void update_particle_list( int level ) { 
	int i, k;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	double pos[nDim];
	int num_parts_to_send[MAX_PROCS];
	int particle_list_to_send[MAX_PROCS];
	int *particle_array_to_send[MAX_PROCS];
	int ipart2, new_cell;
	int proc;
	int collect_level;
	int sfc;

	start_time( UPDATE_PARTS_TIMER );
	start_time( WORK_TIMER );

	/* now move particles from one cell list to another */
	for ( i = 0; i < num_procs; i++ ) {
		num_parts_to_send[i] = 0;
		particle_list_to_send[i] = NULL_PARTICLE;
	}

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( k = 0; k < num_level_cells; k++ ) {
		iter_cell = level_cells[k];

		ipart = cell_particle_list[iter_cell];
		while ( ipart != NULL_PARTICLE ) {
			ipart2 = particle_list_next[ipart];

			sfc = sfc_index_position( particle_x[ipart] );
			proc = processor_owner(sfc);

			if ( proc == local_proc_id ) {
				new_cell = cell_find_position_sfc( sfc, particle_x[ipart] );
				if ( new_cell != iter_cell ) {
					cart_assert( cell_is_local(new_cell) );
					delete_particle( iter_cell, ipart );
					insert_particle( new_cell, ipart );
				}
			} else if ( proc == -1 ) {
				cart_error( "Unable to locate processor for particle %d!", particle_id[ipart]);
			} else {
				delete_particle( iter_cell, ipart );
				particle_list_next[ipart] = particle_list_to_send[proc];
				particle_list_to_send[proc] = ipart;
				num_parts_to_send[proc]++;
			}

			ipart = ipart2;
		}
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );

	start_time( COMMUNICATION_TIMER );
	start_time( UPDATE_PARTS_COMMUNICATION_TIMER );

	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_parts_to_send[proc] > 0 ) {
			particle_array_to_send[proc] = cart_alloc(int, num_parts_to_send[proc]);
			num_parts_to_send[proc] = 0;

			/* add particles that ended up in processor linked list */
			ipart = particle_list_to_send[proc];
			while ( ipart != NULL_PARTICLE ) {
				particle_array_to_send[proc][ num_parts_to_send[proc]++ ] = ipart;
				ipart = particle_list_next[ipart];
			}
		}
	}

	trade_particle_lists( num_parts_to_send, particle_array_to_send, level, FREE_PARTICLE_LISTS );

	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_parts_to_send[proc] > 0 ) {
			cart_free( particle_array_to_send[proc] );
		}
    }

	end_time( UPDATE_PARTS_COMMUNICATION_TIMER );
	end_time( COMMUNICATION_TIMER );
	end_time( UPDATE_PARTS_TIMER );
}

void trade_particle_lists( int num_parts_to_send[MAX_PROCS], int *particle_list_to_send[MAX_PROCS], int trade_level, int free_particle_flag ) {
	int i, j;
	int id_count, part_count;
	int proc, icell, ipart;
	int page_size, parts_page_size;
	int proc_pages_sent;
	int num_pages_sent;
	int num_pages_to_send;
	int num_send_requests;
	int num_requests;
	int num_parts_to_recv[MAX_PROCS];
	int page_count[MAX_PROCS];
	int *send_id;
	int *recv_id[MAX_PROCS];
	double *send_parts;
	double *recv_parts[MAX_PROCS];
	int num_pages_received;

#ifdef GRAVITY
	float *send_potential;
	float *recv_potential[MAX_PROCS];
	MPI_Request recv_pot_requests[MAX_PROCS];
#endif /* GRAVITY */

#ifdef STAR_FORMATION
	int star_page_start;
	int star_vars_sent;
	int star_vars_recv;
	int total_star_vars;
	int star_page_size;
	MPI_Request recv_stars_requests[MAX_PROCS];
	float *send_stars;
	float *recv_stars[MAX_PROCS];

#ifdef STAR_PARTICLE_TYPES
	MPI_Request recv_star_type_requests[MAX_PROCS];
	int *send_star_types;
	int *recv_star_types[MAX_PROCS];
#endif /* STAR_PARTICLE_TYPES */

#ifdef ENRICHMENT
#ifdef ENRICHMENT_SNIa
	#define num_star_vars	5
#else
	#define num_star_vars	4
#endif /* ENRICHMENT_SNIa */
#else
	#define num_star_vars	3
#endif /* ENRICHMENT */
#endif /* STAR_FORMATION */

	#define num_particle_vars	(2+2*nDim)	/* t, dt, x, v */

	int num_request_types;
	MPI_Request *send_requests;
	MPI_Request recv_id_requests[MAX_PROCS];
	MPI_Request recv_parts_requests[MAX_PROCS];
	MPI_Request send_count_requests[MAX_PROCS];
	MPI_Status status;

	cart_assert( trade_level == -1 || ( trade_level >= min_level && trade_level <= max_level ) );

	start_time( TRADE_PARTICLE_TIMER );

	/* use same page size as for I/O, could easily change to different parameter,
	 * doesn't really matter as long as page_size is small relative to memory,
	 * but typical of numbers of particles moved */
	page_size = MIN(65536/num_procs, 1024);
	parts_page_size = num_particle_vars*page_size;

#ifdef STAR_FORMATION
	star_page_size = num_star_vars*page_size;
#endif /* STAR_FORMATION */

	/* compute number of pages we'll need to send */
	num_pages_to_send = 0;
	num_requests = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		num_parts_to_recv[proc] = 0;

		if ( ( trade_level == -1 && proc != local_proc_id ) ||
				( trade_level != -1 &&
				( num_remote_buffers[trade_level][proc] > 0
				|| num_local_buffers[trade_level][proc] > 0 ) ) ) {

			if ( num_parts_to_send[proc] > 0 ) {
				num_pages_to_send += ( num_parts_to_send[proc] - 1 ) / page_size + 1;
			}

			MPI_Irecv( &num_parts_to_recv[proc], 1, MPI_INT, proc, 0,
				mpi.comm.run, &recv_id_requests[num_requests] );
			MPI_Isend( &num_parts_to_send[proc], 1, MPI_INT, proc, 0,
				mpi.comm.run, &send_count_requests[num_requests] );

			num_requests++;
		}
	}

	start_time( TRADE_PARTICLE_COMMUNICATION_TIMER );
	MPI_Waitall( num_requests, send_count_requests, MPI_STATUSES_IGNORE );
	MPI_Waitall( num_requests, recv_id_requests, MPI_STATUSES_IGNORE );
	end_time( TRADE_PARTICLE_COMMUNICATION_TIMER );

	/* set up receives */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_parts_to_recv[proc] > 0 ) {
			recv_id[proc] = cart_alloc(int, page_size );
			recv_parts[proc] = cart_alloc(double, parts_page_size );

			MPI_Irecv( recv_id[proc], page_size, MPI_INT, proc, 0, 
				mpi.comm.run, &recv_id_requests[proc] );
			MPI_Irecv( recv_parts[proc], parts_page_size, MPI_DOUBLE, 
				proc, 0, mpi.comm.run, &recv_parts_requests[proc] );

#ifdef GRAVITY
			recv_potential[proc] = cart_alloc(float, page_size );
			MPI_Irecv( recv_potential[proc], page_size, MPI_FLOAT, proc, 1,
				mpi.comm.run, &recv_pot_requests[proc] );
#endif /* GRAVITY */

#ifdef STAR_FORMATION
			recv_stars[proc] = cart_alloc(float, star_page_size );
			MPI_Irecv( recv_stars[proc], star_page_size, MPI_FLOAT, proc, 0, 
				mpi.comm.run, &recv_stars_requests[proc] );
#ifdef STAR_PARTICLE_TYPES
			recv_star_types[proc] = cart_alloc(int, page_size);
			MPI_Irecv( recv_star_types[proc], page_size, MPI_INT, proc, 1,
				mpi.comm.run, &recv_star_type_requests[proc] );
#endif /* STAR_PARTICLE_TYPES */
#endif /* STAR_FORMATION */

			page_count[proc] = 0;
		} else {
			recv_id_requests[proc] = MPI_REQUEST_NULL;
			recv_parts_requests[proc] = MPI_REQUEST_NULL;

#ifdef GRAVITY
			recv_pot_requests[proc] = MPI_REQUEST_NULL;
#endif /* GRAVITY */
#ifdef STAR_FORMATION
			recv_stars_requests[proc] = MPI_REQUEST_NULL;

#ifdef STAR_PARTICLE_TYPES
			recv_star_type_requests[proc] = MPI_REQUEST_NULL;
#endif /* STAR_PARTICLE_TYPES */
#endif /* STAR_FORMATION */
		}
	}

	/* allocate space for the pages */
	send_id = cart_alloc(int, num_pages_to_send * page_size );
	send_parts = cart_alloc(double, num_pages_to_send * parts_page_size );
	num_request_types = 2;

#ifdef GRAVITY
	send_potential = cart_alloc(float, num_pages_to_send*page_size );
	num_request_types++;
#endif /* GRAVITY */
#ifdef STAR_FORMATION
	/* allocating enough space for all particles to be stars */
	send_stars = cart_alloc(float, num_pages_to_send * star_page_size );
	star_page_start = 0;
	star_vars_sent = 0;
	num_request_types++;

#ifdef STAR_PARTICLE_TYPES
	num_request_types++;
	send_star_types = cart_alloc(int, num_pages_to_send*page_size);
#endif /* STAR_PARTICLE_TYPES */
#endif /* STAR_FORMATION */

	send_requests = cart_alloc(MPI_Request, num_request_types*num_pages_to_send);

	num_pages_sent = 0;
	num_send_requests = 0;

	for ( proc = 0; proc < num_procs; proc++ ) {
		id_count = num_pages_sent*page_size;
		part_count = num_pages_sent*parts_page_size;
		proc_pages_sent = 0;

		for ( i = 0; i < num_parts_to_send[proc]; i++ ) {
			ipart = particle_list_to_send[proc][i];
			cart_assert( ipart == NULL_PARTICLE || ( ipart >= 0 && ipart < num_particles ) );
	
			send_id[id_count] = particle_id[ipart];
	
			send_parts[part_count++] = particle_t[ipart];
			send_parts[part_count++] = particle_dt[ipart];
	
			for ( j = 0; j < nDim; j++ ) {
				send_parts[part_count++] = particle_x[ipart][j];
			}
	
			for ( j = 0; j < nDim; j++ ) {
				send_parts[part_count++] = particle_v[ipart][j];
			}

#ifdef GRAVITY
			send_potential[id_count] = particle_pot[ipart];
#endif /* GRAVITY */

#ifdef STAR_FORMATION
			if ( particle_is_star( ipart ) ) {
				/* pack in star variables */
				send_stars[star_vars_sent++] = particle_mass[ipart];
				send_stars[star_vars_sent++] = star_initial_mass[ipart];
				send_stars[star_vars_sent++] = star_tbirth[ipart];

#ifdef ENRICHMENT
				send_stars[star_vars_sent++] = star_metallicity_II[ipart];
#ifdef ENRICHMENT_SNIa
				send_stars[star_vars_sent++] = star_metallicity_Ia[ipart];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef STAR_PARTICLE_TYPES
				send_star_types[id_count] = star_particle_type[ipart];
#endif /* STAR_PARTICLE_TYPES */
			}
#endif /* STAR_FORMATION */

			id_count++;

			if ( id_count % page_size == 0 || i == num_parts_to_send[proc]-1 ) {
				MPI_Isend( &send_id[num_pages_sent*page_size], 
					id_count - num_pages_sent*page_size, 
					MPI_INT, proc, 2*proc_pages_sent, mpi.comm.run, 
					&send_requests[num_send_requests++] );

				MPI_Isend( &send_parts[num_pages_sent*parts_page_size], 
						part_count - num_pages_sent*parts_page_size, 
						MPI_DOUBLE, proc, proc_pages_sent, mpi.comm.run, 
						&send_requests[num_send_requests++] );

#ifdef GRAVITY
				MPI_Isend( &send_potential[num_pages_sent*page_size],
					id_count - num_pages_sent*page_size,
					MPI_FLOAT, proc, 2*proc_pages_sent+1, mpi.comm.run,
					&send_requests[num_send_requests++] );
#endif /* GRAVITY */

#ifdef STAR_FORMATION
				MPI_Isend( &send_stars[star_page_start], 
						star_vars_sent-star_page_start,
						MPI_FLOAT, proc, 2*proc_pages_sent, mpi.comm.run,
						&send_requests[num_send_requests++] );

#ifdef STAR_PARTICLE_TYPES
				/* sending one int per particle, not just star particles */
				MPI_Isend( &send_star_types[num_pages_sent*page_size], 
						id_count - num_pages_sent*page_size, MPI_INT, proc, 
						2*proc_pages_sent+1, mpi.comm.run, 
						&send_requests[num_send_requests++] ); 
#endif /* STAR_PARTICLE_TYPES */

				star_page_start = star_vars_sent;
#endif /* STAR_FORMATION */

				proc_pages_sent++;
				num_pages_sent++;
			}
		}

		if ( free_particle_flag ) {
			for ( i = 0; i < num_parts_to_send[proc]; i++ ) {
				particle_free( particle_list_to_send[proc][i] );
			}
		}	
	}

	cart_assert( num_pages_sent == num_pages_to_send );

	/* wait for receives and process particles */
	num_pages_received = 0;
	do {
		start_time( TRADE_PARTICLE_COMMUNICATION_TIMER );
		MPI_Waitany( num_procs, recv_id_requests, &proc, &status );
		end_time( TRADE_PARTICLE_COMMUNICATION_TIMER );

		if ( proc != MPI_UNDEFINED ) {
			num_pages_received++;
			MPI_Get_count( &status, MPI_INT, &id_count );

			MPI_Wait( &recv_parts_requests[proc], MPI_STATUS_IGNORE );

#ifdef GRAVITY 
			MPI_Wait( &recv_pot_requests[proc], MPI_STATUS_IGNORE );
#endif /* GRAVITY */

#ifdef STAR_FORMATION
			MPI_Wait( &recv_stars_requests[proc], &status );
			MPI_Get_count( &status, MPI_FLOAT, &total_star_vars );

#ifdef STAR_PARTICLE_TYPES
			MPI_Wait( &recv_star_type_requests[proc], MPI_STATUS_IGNORE );
#endif /* STAR_PARTICLE_TYPES */

			star_vars_recv = 0;
#endif /* STAR_FORMATION */

			/* process received page */
			part_count = 0;
			for ( i = 0; i < id_count; i++ ) {
				ipart = particle_alloc( recv_id[proc][i] );
				cart_assert( ipart >= 0 && ipart < num_particles );
				cart_assert( particle_level[ipart] == FREE_PARTICLE_LEVEL );

				particle_t[ipart] = recv_parts[proc][part_count++];
				particle_dt[ipart] = recv_parts[proc][part_count++];

				for ( j = 0; j < nDim; j++ ) {
					particle_x[ipart][j] = recv_parts[proc][part_count++];
				}

				for ( j = 0; j < nDim; j++ ) {
					particle_v[ipart][j] = recv_parts[proc][part_count++];
				}

#ifdef GRAVITY 
				particle_pot[ipart] = recv_potential[proc][i];
#endif /* GRAVITY */

#ifdef STAR_FORMATION
				if ( particle_id_is_star( recv_id[proc][i] ) ) {
					cart_assert( ipart >= 0 && ipart < num_star_particles );
					cart_assert( particle_is_star(ipart) );

					/* unpack star variables */
					particle_mass[ipart] = recv_stars[proc][star_vars_recv++];
					star_initial_mass[ipart] = recv_stars[proc][star_vars_recv++];
					star_tbirth[ipart] = recv_stars[proc][star_vars_recv++];

#ifdef ENRICHMENT
					star_metallicity_II[ipart] = recv_stars[proc][star_vars_recv++];
#ifdef ENRICHMENT_SNIa
					star_metallicity_Ia[ipart] = recv_stars[proc][star_vars_recv++];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

#ifdef STAR_PARTICLE_TYPES
					star_particle_type[ipart] = recv_star_types[proc][i];
#endif /* STAR_PARTICLE_TYPES */
				
				} else {
					particle_mass[ipart] = particle_species_mass[ particle_species( recv_id[proc][i] ) ];
				}
#else
				particle_mass[ipart] = particle_species_mass[ particle_species( recv_id[proc][i] ) ];
#endif /* STAR_FORMATION */

				icell = cell_find_position( particle_x[ipart] );
				cart_assert( icell != -1 );
				insert_particle( icell, ipart );
			}

#ifdef STAR_FORMATION
			/* ensure we unpacked all stars */
			cart_assert( star_vars_recv == total_star_vars );
#endif /* STAR_FORMATION */

			num_parts_to_recv[proc] -= id_count;

			/* if we received a full page, set up to receive a new one */
			if ( num_parts_to_recv[proc] > 0 ) {
				page_count[proc]++;
				MPI_Irecv( recv_id[proc], page_size, MPI_INT, proc, 2*page_count[proc], 
						mpi.comm.run, &recv_id_requests[proc] );
				MPI_Irecv( recv_parts[proc], parts_page_size, MPI_DOUBLE, proc, 
						page_count[proc], mpi.comm.run, &recv_parts_requests[proc] );

#ifdef GRAVITY
				MPI_Irecv( recv_potential[proc], page_size, MPI_FLOAT, proc,
						2*page_count[proc]+1, mpi.comm.run, &recv_pot_requests[proc] );
#endif

#ifdef STAR_FORMATION
				MPI_Irecv( recv_stars[proc], star_page_size, MPI_FLOAT, proc, 
						2*page_count[proc], mpi.comm.run, &recv_stars_requests[proc] );

#ifdef STAR_PARTICLE_TYPES
				MPI_Irecv( recv_star_types[proc], page_size, MPI_INT, proc, 2*page_count[proc]+1,
						mpi.comm.run, &recv_star_type_requests[proc] );
#endif /* STAR_PARTICLE_TYPES */
#endif /* STAR_FORMATION */
			} else {
				cart_free( recv_id[proc] );
				cart_free( recv_parts[proc] );
#ifdef GRAVITY
				cart_free( recv_potential[proc] );
#endif /* GRAVITY */
#ifdef STAR_FORMATION
				cart_free( recv_stars[proc] );

#ifdef STAR_PARTICLE_TYPES
				cart_free( recv_star_types[proc] );
#endif /* STAR_PARTICLE_TYPES */
#endif /* STAR_FORMATION */
			}
		}
	} while ( proc != MPI_UNDEFINED );

	/* wait for sends to complete */
	start_time( TRADE_PARTICLE_COMMUNICATION_TIMER );
	MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	end_time( TRADE_PARTICLE_COMMUNICATION_TIMER );

	/* de-allocate send buffers */
	cart_free( send_id );
	cart_free( send_parts );
	cart_free( send_requests );
#ifdef GRAVITY
	cart_free( send_potential );
#endif /* GRAVITY */
#ifdef STAR_FORMATION
	cart_free( send_stars );
#ifdef STAR_PARTICLE_TYPES
	cart_free( send_star_types );
#endif /* STAR_PARTICLE_TYPES */
#endif /* STAR_FORMATION */

	end_time( TRADE_PARTICLE_TIMER );
}

void rebuild_particle_list() 
/* recreates particle_list_prev links based on current particle_list_next lists */
{                                                                                                       
	int i;
	int ipart;                                                                                                                                      
	int level, num_level_cells;
	int *level_cells;
	int icell;
	int prev;

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			prev = NULL_PARTICLE;
			ipart = cell_particle_list[icell];
			while ( ipart != NULL_PARTICLE ) {
				particle_list_prev[ipart] = prev;
				prev = ipart;
				ipart = particle_list_next[ipart];
			}
		}

		cart_free( level_cells );
	}
}

void build_particle_list() {
	int i, j;
	int icell;

	cart_debug("build_particle_list()");

	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_id[i] != NULL_PARTICLE ) {
			/* enforce periodicity */
			for ( j = 0; j < nDim; j++ ) {
				if ( particle_x[i][j] >= (double)num_grid ) {
					particle_x[i][j] -= (double)num_grid;
				}

				if ( particle_x[i][j] < 0.0 ) {
					particle_x[i][j] += (double)num_grid;
				}
			}

			/* find cell particle belongs to */
			icell = cell_find_position( particle_x[i] );

			/* some consistency checks */
			cart_assert( icell >= 0 && icell < num_cells );
			cart_assert( cell_is_leaf( icell ) );
			if ( !cell_contains_position(icell, particle_x[i]) ) {
				cart_debug("%d %e %e %e", icell, particle_x[i][0], particle_x[i][1], particle_x[i][2] );
				cart_error("Error in building the particle list");
			}

			/* insert particle into cell list */
			insert_particle( icell, i );
		}
	}

	particle_list_enabled = 1;
}

int particle_alloc( int id ) { 
	int ipart;
	int i;

#ifdef STAR_FORMATION
	if ( particle_id_is_star(id) ) {
		if ( free_star_particle_list == NULL_PARTICLE ) {
			/* search for first free star particle */
			while ( next_free_star_particle < num_star_particles &&
					particle_level[next_free_star_particle] != FREE_PARTICLE_LEVEL ) {
				next_free_star_particle++;
			}

			if ( next_free_star_particle < num_star_particles ) {
				ipart = next_free_star_particle;
				next_free_star_particle++;
			} else {
				/* block from 0->num_star_particles completely filled, try to find a non-star */
				for ( i = 0; i < num_star_particles; i++ ) {
					if ( !particle_is_star(i) ) {
						/* move non-star particle to make room for new star */
						ipart = particle_alloc( particle_id[i] );
						num_local_particles--;
						particle_move( i, ipart );
						ipart = i;
						break;
					}
				}

				if ( i == num_star_particles ) {
                                        cart_error("Ran out of star particles %d, increase num_star_particles!", i);
				}
			}
		} else {
			/* take off the recently removed stack */
			ipart = free_star_particle_list;
			free_star_particle_list = particle_list_next[free_star_particle_list];
		}

		cart_assert( ipart >= 0 && ipart < num_star_particles );
		num_local_star_particles++;
	} else {
		if ( free_particle_list == NULL_PARTICLE ) {
			if ( next_free_particle >= num_particles ) {
				/* ran out of normal particles, start using stars */
				if ( free_star_particle_list == NULL_PARTICLE ) {
					while ( next_free_star_particle < num_star_particles &&
							particle_level[next_free_star_particle] != FREE_PARTICLE_LEVEL ) {
						next_free_star_particle++;
					}

					if ( next_free_star_particle >= num_star_particles ) {
                                                cart_error("Ran out of particles (next=%d), increase num_star_particles!",next_free_star_particle );
					} else {
						ipart = next_free_star_particle;
						next_free_star_particle++;
					}
				} else {
					ipart = free_star_particle_list;
					free_star_particle_list = particle_list_next[free_star_particle_list];
				}
			} else {
				ipart = next_free_particle;
				next_free_particle++;
			}
		} else {
			ipart = free_particle_list;
			free_particle_list = particle_list_next[free_particle_list];
		}
	}
#else 
	if ( free_particle_list == NULL_PARTICLE ) {
		if ( num_local_particles >= num_particles ) {
			/* generate an error, ran out of particles */
                        cart_error("Ran out of local particles %d, increase num_particles!", num_local_particles);
		}

		ipart = next_free_particle;
		next_free_particle++;
	} else {
		ipart = free_particle_list;
		free_particle_list = particle_list_next[free_particle_list];
	}
#endif /* STAR_FORMATION */

	cart_assert( ipart >= 0 && ipart < num_particles );
	cart_assert( particle_level[ipart] == FREE_PARTICLE_LEVEL );
	cart_assert( particle_id[ipart] == NULL_PARTICLE );

	particle_id[ipart] = id;
	num_local_particles++;

	return ipart;
}

void particle_move( int ipart_old, int ipart_new ) {
	int icell;
	int i;

	cart_assert( particle_level[ipart_old] != FREE_PARTICLE_LEVEL );
	cart_assert( particle_level[ipart_new] == FREE_PARTICLE_LEVEL );

	particle_id[ipart_new] = particle_id[ipart_old];
	particle_level[ipart_new] = particle_level[ipart_old];
	particle_t[ipart_new] = particle_t[ipart_old];
	particle_dt[ipart_new] = particle_dt[ipart_old];
	particle_mass[ipart_new] = particle_mass[ipart_old];

#ifdef GRAVITY
	particle_pot[ipart_new] = particle_pot[ipart_old];
#endif /* GRAVITY */

	for ( i = 0; i < nDim; i++ ) {
		particle_x[ipart_new][i] = particle_x[ipart_old][i];
	}

	for ( i = 0; i < nDim; i++ ) {
		particle_v[ipart_new][i] = particle_x[ipart_old][i];
	}

	if ( particle_list_enabled ) {
		if ( particle_list_next[ipart_old] != NULL_PARTICLE ) {
			cart_assert( particle_list_prev[ particle_list_next[ipart_old] ] == ipart_old );
			particle_list_prev[ particle_list_next[ipart_old] ] = ipart_new;
		}

		if ( particle_list_prev[ipart_old] == NULL_PARTICLE ) {
			/* this should be first particle in cell list */
			icell = cell_find_position( particle_x[ipart_old] );

			cart_assert( cell_particle_list[icell] == ipart_old );
			cell_particle_list[icell] = ipart_new;
		} else {
			cart_assert( particle_list_next[ particle_list_prev[ipart_old] ] == ipart_old );
			particle_list_next[ particle_list_prev[ipart_old] ] = ipart_new;
		}

		particle_list_next[ipart_new] = particle_list_next[ipart_old];
		particle_list_prev[ipart_new] = particle_list_prev[ipart_old];
	}

#ifdef STAR_FORMATION
	/* move star variables here */
	if ( particle_is_star(ipart_old) ) {
		cart_assert( ipart_new < num_star_particles );

		star_tbirth[ipart_new] = star_tbirth[ipart_old];
		star_initial_mass[ipart_new] = star_initial_mass[ipart_old];
		
#ifdef ENRICHMENT
		star_metallicity_II[ipart_new] = star_metallicity_II[ipart_old];
#ifdef ENRICHMENT_SNIa
		star_metallicity_Ia[ipart_new] = star_metallicity_Ia[ipart_old];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef STAR_PARTICLE_TYPES
        star_particle_type[ipart_new] = star_particle_type[ipart_old];
        star_particle_type[ipart_old] = STAR_TYPE_DELETED;
#endif /* STAR_PARTICLE_TYPES */	
	}
#endif /* STAR_FORMATION */

	particle_id[ipart_old] = NULL_PARTICLE;
	particle_level[ipart_old] = FREE_PARTICLE_LEVEL;
}

void particle_free( int ipart ) {
	cart_assert( ipart >= 0 && ipart < num_particles );
	cart_assert( particle_level[ipart] != FREE_PARTICLE_LEVEL );

#ifdef STAR_FORMATION
	if ( ipart < num_star_particles ) {
		particle_list_next[ipart] = free_star_particle_list;
		free_star_particle_list = ipart;

		if ( particle_is_star(ipart) ) {
#ifdef STAR_PARTICLE_TYPES
            star_particle_type[ipart] = STAR_TYPE_DELETED;
#endif /* STAR_PARTICLE_TYPES */

			num_local_star_particles--;
		}
	} else {
		/* not a star, add to normal list */
		particle_list_next[ipart] = free_particle_list;
		free_particle_list = ipart;
	}
#else
	/* not a star, add to normal list */
	particle_list_next[ipart] = free_particle_list;
	free_particle_list = ipart;
#endif /* STAR_FORMATION */

	particle_level[ipart] = FREE_PARTICLE_LEVEL;
	particle_id[ipart] = NULL_PARTICLE;

	num_local_particles--;
}

void particle_list_free( int ihead ) {
	int last, next;

	cart_assert( ihead != NULL_PARTICLE );
	last = ihead;

	while ( last != NULL_PARTICLE ) {
		next = particle_list_next[last];

#ifdef STAR_FORMATION
		if ( last < num_star_particles ) {
			particle_list_next[last] = free_star_particle_list;
			free_star_particle_list = last;

			if ( particle_is_star(last) ) {
				num_local_star_particles--;
			}
		} else {
			particle_list_next[last] = free_particle_list;
			free_particle_list = last;
		}
#else
		particle_list_next[last] = free_particle_list;
		free_particle_list = last;
#endif /* STAR_FORMATION */

		particle_level[last] = FREE_PARTICLE_LEVEL;
		particle_id[last] = NULL_PARTICLE;
		num_local_particles--;

		last = next;
	}
}

void split_particle_list( int cell ) {
	int i;
	int part;
	int next;
	int child;
	double pos[nDim];

	cart_assert( cell >= 0 && cell < num_cells );
	cart_assert( cell_is_refined(cell) );

	part = cell_particle_list[cell];
	cell_particle_list[cell] = NULL_PARTICLE;

	cell_center_position( cell, pos );

	while ( part != NULL_PARTICLE ) {
		next = particle_list_next[part];

		/* find which child contains this particle */
		child = 0;
		for ( i = 0; i < nDim; i++ ) {
			if ( particle_x[part][i] >= pos[i] ) {
				child += (1<<i);
			}
		}

		insert_particle( cell_child( cell, child ), part );
		part = next;
	}
}

void join_particle_list( int cell ) {
	int i;
	int level;
	int part;
	int head;

	cart_assert( cell >= 0 && cell < num_cells );
	cart_assert( cell_is_refined(cell) );
	
	level = cell_level(cell);

	/* find children with particles */
	part = NULL_PARTICLE;
	for ( i = 0; i < num_children; i++ ) {
		head = cell_particle_list[ cell_child( cell, i ) ];
		cell_particle_list[ cell_child( cell, i ) ] = NULL_PARTICLE;

		if ( head != NULL_PARTICLE ) {
			if ( part == NULL_PARTICLE ) {
				cell_particle_list[cell] = head;
			} else {
				particle_list_next[part] = head;
				particle_list_prev[head] = part;
			}

			part = head;
			particle_level[part] = level;
			while ( particle_list_next[part] != NULL_PARTICLE ) {
				part = particle_list_next[part];
				particle_level[part] = level;
			}
		}
	}
}

void insert_particle( int cell, int part ) {
	int head;

	cart_assert( cell >= 0 && cell < num_cells );
	cart_assert( cell_is_leaf( cell ) );
	cart_assert( part >= 0 && part < num_particles );
	cart_assert( cell_level(cell) >= min_level && cell_level(cell) <= max_level );

	head = cell_particle_list[cell];

	particle_list_prev[part] = NULL_PARTICLE;
	particle_list_next[part] = head;

	if ( head != NULL_PARTICLE ) {
		cart_assert( particle_list_prev[head] == NULL_PARTICLE );
		particle_list_prev[head] = part;
	}
	
	cell_particle_list[cell] = part;
	particle_level[part] = cell_level(cell);
}

void delete_particle( int icell, int part ) {
	int next, prev;

	cart_assert( part >= 0 && part < num_particles );
	cart_assert( particle_level[part] != FREE_PARTICLE_LEVEL );

	prev = particle_list_prev[part];
	next = particle_list_next[part];

	if ( prev == NULL_PARTICLE ) {
		cart_assert( cell_particle_list[icell] == part );
		cell_particle_list[icell] = next;
	} else {
		particle_list_next[prev] = next;
	}

	if ( next != NULL_PARTICLE ) {
		particle_list_prev[next] = prev;
	}
}

int particle_species( int id ) {
	int specie = 0;

	for ( specie = 0; specie < num_particle_species; specie++ ) {
		if ( id < particle_species_indices[specie+1] ) {
			break;
		}
	}

	return specie;
}

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
#ifdef REFINEMENT
void get_refinement_region(){
	int i,j;
	float refmin[nDim];
	float refmax[nDim];

	/* Doug (11/29/2009): necessary to properly set particle timestep 
	 * (not certain abox and auni are required here) */
	/* NG: can NOT be used here, hydro vars may not have been read yet! */
	//set_timestepping_scheme();
	
	if ( spatially_limited_refinement ) {
		for ( i = 0; i < nDim; i++ ) {
			refmin[i] = refinement_volume_min[i];
			refmax[i] = refinement_volume_max[i];
		}

		for ( j = 0; j < num_particles; j++ ) {
			if ( particle_level[j] != FREE_PARTICLE_LEVEL &&
						particle_id[j] < particle_species_indices[1] ) {
				for ( i = 0; i < nDim; i++ ) {
					if ( particle_x[j][i] < refmin[i] ) {
						refmin[i] = particle_x[j][i];
					}

					if ( particle_x[j][i] > refmax[i] ) {
						refmax[i] = particle_x[j][i];
					}
				}
			}
		}

		for ( i = 0; i < nDim; i++ ) {
			refmin[i] = floor(refmin[i]);
			refmax[i] = ceil(refmax[i]);
		}

		MPI_Allreduce( refmin, refinement_volume_min, nDim, MPI_FLOAT, MPI_MIN, mpi.comm.run );
		MPI_Allreduce( refmax, refinement_volume_max, nDim, MPI_FLOAT, MPI_MAX, mpi.comm.run );

		for ( i = 0; i < nDim; i++ ) {
			cart_debug("refinement_volume[%u] = %e %e", i, refinement_volume_min[i], refinement_volume_max[i] );
		}
	}
}

void build_refinement_region(int do_load_balance){
	int j;
	int level, cell;
	int total_cells_per_level[max_level-min_level+1];

	/* do initial refinement */
	level = min_level;
	total_cells_per_level[min_level] = num_root_cells;
	while ( level < max_level && total_cells_per_level[level] > 0 ) {
		cart_debug("assigning density to level %u", level );
		assign_density(level);
		cart_debug("refining level %u, num_cells_per_level = %d", level, num_cells_per_level[level] );
		modify( level, OP_REFINE );
		cart_debug("done refining level %u, created %u new cells", 
				level, num_cells_per_level[level+1] );
		MPI_Allreduce( &num_cells_per_level[level+1], &total_cells_per_level[level+1], 1, MPI_INT, MPI_SUM, mpi.comm.run );
		level++;
	
		if ( local_proc_id == MASTER_NODE ) {
			cart_debug("level %u: %u cells", level, total_cells_per_level[level] );
		}

		if(total_cells_per_level[level] > 0) {
			for(j=0; j<num_particles; j++) { 
				if ( particle_level[j] != FREE_PARTICLE_LEVEL) {
					cell = cell_find_position_above_level(level,particle_x[j]);
					particle_level[j] = cell_level(cell);
				}
			}
			if(do_load_balance){
				load_balance();
			}
		}
	}
}
#endif /* REFINEMENT */

void build_mesh() {
	/* Doug (11/29/2009): necessary to properly set particle timestep 
	 * (not certain abox and auni are required here) */
	/* NG: can NOT be used here, hydro vars may not have been read yet! */
	//set_timestepping_scheme();

#ifdef REFINEMENT
	get_refinement_region();
    build_cell_buffer();
    repair_neighbors();
	build_refinement_region(1);
#else
	build_cell_buffer();
	repair_neighbors();
	load_balance();
#endif /* REFINEMENT */
}

#endif /* GRAVITY || RADIATIVE_TRANSFER */

#endif /* PARTICLES */



