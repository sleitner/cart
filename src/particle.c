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
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "sfc.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


int num_row;

double particle_t[num_particles];
double particle_dt[num_particles];
double particle_x[num_particles][nDim];
double particle_v[num_particles][nDim];

#ifdef GRAVITY
float particle_pot[num_particles];
#endif /* GRAVITY */

int particle_level[num_particles];
float particle_mass[num_particles];
int particle_id[num_particles];
int particle_list_next[num_particles];
int particle_list_prev[num_particles];

/* particle species */
int num_particle_species;
float particle_species_mass[MAX_PARTICLE_SPECIES];
int particle_species_num[MAX_PARTICLE_SPECIES];
int particle_species_indices[MAX_PARTICLE_SPECIES+1];

/* variables for logging energy */
double tintg;
double ekin;
double ekin1;
double ekin2;
double au0;
double aeu0;
double ap0, ap1;

int cell_particle_list[num_cells];

int num_local_particles = 0;
long num_particles_total = 0;

int next_free_particle = 0;
int free_particle_list = NULL_PARTICLE;

#ifdef STARFORM
int next_free_star_particle = 0;
int free_star_particle_list = NULL_PARTICLE;
#endif /* STARFORM */

int particle_list_enabled = 0;

void init_particles() { 
	int i;

#pragma omp parallel for default(none), private(i), shared(cell_particle_list)
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

#ifdef STARFORM
	free_star_particle_list = NULL_PARTICLE;
	next_free_star_particle = 0;
	next_free_particle = num_star_particles;
	num_local_star_particles = 0;
#else 
	next_free_particle = 0;
#endif /* STARFORM */

	particle_list_enabled = 0;
}


void move_particles( int level ) {
	int i, m;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	int icell, icell_orig;
	int level1;
	int child;
	double pos[nDim];
	int found;
	double t_next;
	int c[num_children];
	double diff1, diff2, diff3;
	double pt3, pd3;
	double t1,t2,t3,d1,d2,d3;
	double t2t1, t2d1, d2t1, d2d1;
	double x, y, z, vx, vy, vz, ax, ay, az;
	double t3t2t1, t3t2d1, t3d2t1, t3d2d1;
	double d3t2t1, d3t2d1, d3d2t1, d3d2d1;
	double pconst;
	double delta_t;

	cart_assert( level >= min_level && level <= max_level );
        
	start_time( MOVE_PARTS_TIMER );
	start_time( WORK_TIMER ); 

#ifdef STARFORM
	setup_star_formation_feedback(level);
#endif /* STARFORM */

	t_next = tl[level] + dtl[level];

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#ifdef GRAVITY
#pragma omp parallel for default(none), private(iter_cell,ipart,x,y,z,icell,level1,found,icell_orig,pos,child,i,c,diff1,diff2,diff3,d1,d2,d3,t1,t2,t3,pconst,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1,d3t2t1,d3t2d1,d3d2t1,d3d2d1,ax,ay,az,vx,vy,vz,delta_t,t2t1,t2d1,d2t1,d2d1), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,particle_x,t_next,particle_dt,particle_v,particle_id,particle_species_indices,num_particle_species,dtl,particle_list_next,cell_size_inverse,tl,particle_pot,particle_mass,cell_vars,neighbor_moves), schedule(dynamic)
#else
#pragma omp parallel for default(none), private(iter_cell,ipart,x,y,z,icell,level1,found,icell_orig,pos,child,i,c,diff1,diff2,diff3,d1,d2,d3,t1,t2,t3,pconst,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1,d3t2t1,d3t2d1,d3d2t1,d3d2d1,ax,ay,az,vx,vy,vz,delta_t), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,particle_x,t_next,particle_dt,particle_v,particle_id,particle_species_indices,num_particle_species,dtl,particle_list_next,neighbor_moves)
#endif
	for ( m = 0; m < num_level_cells; m++ ) {
		iter_cell = level_cells[m];

		ipart = cell_particle_list[iter_cell];
		while ( ipart != NULL_PARTICLE ) {
			cart_assert( ipart >= 0 && ipart < num_particles );
			cart_assert( particle_level[ipart] == level );
		
			if ( particle_t[ipart] < t_next - 0.5*dtl[max_level] ) {
				x = particle_x[ipart][0];
				y = particle_x[ipart][1];
				z = particle_x[ipart][2];

				icell = iter_cell;
				level1 = level;
				do {
					found = 1;
					icell_orig = icell;
					cart_assert( icell != NULL_OCT );

					cell_center_position( icell, pos );

					/* find lower leftmost cell */
					child = 0;
					for ( i = 0; i < nDim; i++ ) {
						if ( particle_x[ipart][i] >= pos[i] ) {
							child += (1<<i);
						}
					}
					cart_assert( child >= 0 && child < num_children );

					for ( i = 0; i < nDim; i++ ) {
						if ( neighbor_moves[child][i] == -1 ) {
							break;
						} else {
							icell = cell_neighbor(icell, neighbor_moves[child][i] );
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
	
						for ( i = 1; i < num_children; i++ ) {
							if ( cell_level(c[i]) != level1 ) {
								icell = cell_parent_cell(icell_orig);
								level1 = level1 - 1;
								cart_assert( icell != NULL_OCT );
								found = 0;
								break;
							}
						}
					}
				} while ( !found );
	
				for ( i = 0; i < num_children; i++ ) {
					cart_assert( cell_level(c[i]) == level1 );
				}

				cart_assert( c[0] != NULL_OCT );
				cell_center_position( c[0], pos );

#ifdef GRAVITY
				/* now we have the level on which this particle will move */
				diff1 = pos[0] - x;
				if ( fabs(diff1) > (double)(num_grid/2) ) {
					if ( diff1 > 0.0 ) {
						diff1 -= (double)(num_grid);
					} else {
						diff1 += (double)(num_grid);
					}
				}
				d1 = fabs(diff1) * cell_size_inverse[level1];
				cart_assert( d1 >= 0.0 && d1 <= 1.0 );

				diff2 = pos[1] - y;
				if ( fabs(diff2) > (double)(num_grid/2) ) {
					if ( diff2 > 0.0 ) {
						diff2 -= (double)(num_grid);
					} else {
						diff2 += (double)(num_grid);
					}
				}
				d2 = fabs(diff2) * cell_size_inverse[level1];

				diff3 = pos[2] - z;
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

				/* a2cs term in HART is computed in accel */
				pconst = tl[level] - particle_t[ipart] + 0.5*( dtl[level] + particle_dt[ipart] );

				pt3 = pconst * t3;
				pd3 = pconst * d3;

				t3t2t1 = pt3 * t2t1;
				t3t2d1 = pt3 * t2d1;
				t3d2t1 = pt3 * d2t1;
				t3d2d1 = pt3 * d2d1;
				d3t2t1 = pd3 * t2t1;
				d3t2d1 = pd3 * t2d1;
				d3d2t1 = pd3 * d2t1;
				d3d2d1 = pd3 * d2d1;

				particle_pot[ipart] = particle_mass[ipart] * 
						(	t3t2t1 * cell_potential(c[0]) +
							t3t2d1 * cell_potential(c[1]) +
							t3d2t1 * cell_potential(c[2]) +
							t3d2d1 * cell_potential(c[3]) +
							d3t2t1 * cell_potential(c[4]) +
							d3t2d1 * cell_potential(c[5]) +
							d3d2t1 * cell_potential(c[6]) +
							d3d2d1 * cell_potential(c[7]) ) / pconst;

				ax = 	t3t2t1 * cell_accel(c[0], 0) +
					t3t2d1 * cell_accel(c[1], 0) +
					t3d2t1 * cell_accel(c[2], 0) +
					t3d2d1 * cell_accel(c[3], 0) + 
					d3t2t1 * cell_accel(c[4], 0) + 
					d3t2d1 * cell_accel(c[5], 0) + 
					d3d2t1 * cell_accel(c[6], 0) + 
					d3d2d1 * cell_accel(c[7], 0);
	
				ay =    t3t2t1 * cell_accel(c[0], 1) +
					t3t2d1 * cell_accel(c[1], 1) +
					t3d2t1 * cell_accel(c[2], 1) +	
					t3d2d1 * cell_accel(c[3], 1) +
					d3t2t1 * cell_accel(c[4], 1) +
					d3t2d1 * cell_accel(c[5], 1) +
					d3d2t1 * cell_accel(c[6], 1) +
					d3d2d1 * cell_accel(c[7], 1);
	
				az =    t3t2t1 * cell_accel(c[0], 2) +
					t3t2d1 * cell_accel(c[1], 2) +
					t3d2t1 * cell_accel(c[2], 2) +
					t3d2d1 * cell_accel(c[3], 2) +
					d3t2t1 * cell_accel(c[4], 2) +
					d3t2d1 * cell_accel(c[5], 2) +
					d3d2t1 * cell_accel(c[6], 2) +
					d3d2d1 * cell_accel(c[7], 2);

				vx = particle_v[ipart][0] + ax;
				vy = particle_v[ipart][1] + ay;
				vz = particle_v[ipart][2] + az;
#else
				vx = particle_v[ipart][0];
				vy = particle_v[ipart][1];
				vz = particle_v[ipart][2];
#endif /* GRAVITY */
	
				delta_t = t_next - particle_t[ipart];

#ifdef STARFORM
				/* do feedback, enrichment, etc */
				if ( particle_is_star(ipart) ) {
				    stellar_feedback(level,iter_cell,ipart,delta_t,t_next,vx,vy,vz);
				}
#endif /* STARFORM */

				x += vx * delta_t;
				y += vy * delta_t;
				z += vz * delta_t;
	
				if ( x < 0.0 ) {
					x += (double)(num_grid);		
				} else if ( x >= (double)(num_grid) ) {
					x -= (double)(num_grid);
				}
	
				if ( y < 0.0 ) {
					y += (double)(num_grid);
				} else if ( y >= (double)(num_grid) ) {
					y -= (double)(num_grid);
				}
	
				if ( z < 0.0 ) {
					z += (double)(num_grid);
				} else if ( z >= (double)(num_grid) ) {
					z -= (double)(num_grid);
				}
	
				particle_x[ipart][0] = x;
				particle_x[ipart][1] = y;
				particle_x[ipart][2] = z;
	
				particle_v[ipart][0] = vx;
				particle_v[ipart][1] = vy;
				particle_v[ipart][2] = vz;
	
				particle_t[ipart] = t_next;
				particle_dt[ipart] = dtl[level];
			}
		
			ipart = particle_list_next[ipart];
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
	end_time( MOVE_PARTS_TIMER );
}

void update_particle_list( int level ) {
	int i, k;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	int coords[nDim];
	int num_parts_to_send[MAX_PROCS];
	int particle_list_to_send[MAX_PROCS];
	int ipart2, new_cell;
	int last_part;
	int proc;
	int collect_level;

	int min_level_modified, max_level_modified;

	start_time( UPDATE_PARTS_TIMER );
	start_time( WORK_TIMER );

	min_level_modified = max_level_modified = level;

	/* now move particles from one cell list to another */
	for ( i = 0; i < num_procs; i++ ) {
		num_parts_to_send[i] = 0;
		particle_list_to_send[i] = NULL_PARTICLE;
	}

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( k = 0; k < num_level_cells; k++ ) {
		iter_cell = level_cells[k];

		if ( cell_is_leaf(iter_cell) ) {
			ipart = cell_particle_list[iter_cell];

			while ( ipart != NULL_PARTICLE ) {
				ipart2 = particle_list_next[ipart];
				new_cell = cell_find_position( particle_x[ipart] );

				if ( new_cell == NULL_OCT ) {
					for ( i = 0; i < nDim; i++ ) {
						coords[i] = (int)(particle_x[ipart][i]);
					}
					proc = processor_owner( sfc_index( coords ) );

					delete_particle( iter_cell, ipart );
					particle_list_next[ipart] = particle_list_to_send[proc];
					particle_list_to_send[proc] = ipart;
				} else { 
					if ( new_cell != iter_cell ) {
						if ( cell_level(new_cell) < min_level_modified ) {
							min_level_modified = cell_level( new_cell );
						}
	
						if ( cell_level(new_cell) > max_level_modified ) {
							max_level_modified = cell_level( new_cell );
						}
	
						delete_particle( iter_cell, ipart );
						insert_particle( new_cell, ipart );
					}
				}

				ipart = ipart2;
			}
		}
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );

	start_time( COMMUNICATION_TIMER );
	start_time( UPDATE_PARTS_COMMUNICATION_TIMER );

	/* now collect particles which ended up in buffer cells */
	for ( collect_level = min_level_modified; collect_level <= max_level_modified; collect_level++ ) {
		select_level( collect_level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			iter_cell = level_cells[i];

			if ( cell_particle_list[iter_cell] != NULL_PARTICLE ) {
				proc = processor_owner( cell_parent_root_sfc( iter_cell ) );

				last_part = cell_particle_list[iter_cell];
				particle_level[last_part] = FREE_PARTICLE_LEVEL;
				num_parts_to_send[proc]++;

				while ( particle_list_next[last_part] != NULL_PARTICLE ) {
					last_part = particle_list_next[last_part];
					particle_level[last_part] = FREE_PARTICLE_LEVEL;
					num_parts_to_send[proc]++;
				}

				particle_list_next[last_part] = particle_list_to_send[proc];
				particle_list_to_send[proc] = cell_particle_list[iter_cell];

				cell_particle_list[iter_cell] = NULL_PARTICLE;
			}
		}
		cart_free( level_cells );
	}

	trade_particle_lists( num_parts_to_send, particle_list_to_send, level );

	end_time( UPDATE_PARTS_COMMUNICATION_TIMER );
	end_time( COMMUNICATION_TIMER );
	end_time( UPDATE_PARTS_TIMER );
}

void trade_particle_lists( int *num_parts_to_send, int *particle_list_to_send, int trade_level ) {
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

#ifdef STARFORM
	int star_page_start;
	int star_vars_sent;
	int star_vars_recv;
	int total_star_vars;
	int star_page_size;
	MPI_Request recv_stars_requests[MAX_PROCS];
	float *send_stars;
	float *recv_stars[MAX_PROCS];

#ifdef ENRICH
#ifdef ENRICH_SNIa
	#define num_star_vars	5
#else
	#define num_star_vars	4
#endif /* ENRICH_SNIa */
#else
	#define num_star_vars	3
#endif /* ENRICH */
#endif /* STARFORM */

	#define num_particle_vars	(2+2*nDim)	/* t, dt, x, v */

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
	page_size = min(num_row*num_row/num_procs, 1024);
	parts_page_size = num_particle_vars*page_size;

#ifdef STARFORM
	star_page_size = num_star_vars*page_size;
#endif /* STARFORM */

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
				MPI_COMM_WORLD, &recv_id_requests[num_requests] );
			MPI_Isend( &num_parts_to_send[proc], 1, MPI_INT, proc, 0,
				MPI_COMM_WORLD, &send_count_requests[num_requests] );

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
				MPI_COMM_WORLD, &recv_id_requests[proc] );
			MPI_Irecv( recv_parts[proc], parts_page_size, MPI_DOUBLE, 
				proc, 0, MPI_COMM_WORLD, &recv_parts_requests[proc] );

#ifdef STARFORM
			recv_stars[proc] = cart_alloc(float, star_page_size );
			MPI_Irecv( recv_stars[proc], star_page_size, MPI_FLOAT, proc, 0, 
				MPI_COMM_WORLD, &recv_stars_requests[proc] );
#endif /* STARFORM */

			page_count[proc] = 0;
		} else {
			recv_id_requests[proc] = MPI_REQUEST_NULL;
			recv_parts_requests[proc] = MPI_REQUEST_NULL;

#ifdef STARFORM
			recv_stars_requests[proc] = MPI_REQUEST_NULL;
#endif /* STARFORM */
		}
	}

	/* allocate space for the pages */
	send_id = cart_alloc(int, num_pages_to_send * page_size );
	send_parts = cart_alloc(double, num_pages_to_send * parts_page_size );

#ifdef STARFORM
	/* allocating enough space for all particles to be stars */
	send_stars = cart_alloc(float, num_pages_to_send * star_page_size );
	send_requests = cart_alloc(MPI_Request, 3*num_pages_to_send );
	star_page_start = 0;
	star_vars_sent = 0;
#else
	send_requests = cart_alloc(MPI_Request, 2*num_pages_to_send );
#endif /* STARFORM */	

	num_pages_sent = 0;
	num_send_requests = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_parts_to_send[proc] > 0 ) {
			part_count = 0;
			id_count = 0;
			proc_pages_sent = 0;

			ipart = particle_list_to_send[proc];
			cart_assert( ipart == NULL_PARTICLE || ( ipart >= 0 && ipart < num_particles ) );

			while ( ipart != NULL_PARTICLE ) {
				send_id[page_size*num_pages_sent+id_count++] = particle_id[ipart];
	
				send_parts[parts_page_size*num_pages_sent+part_count++] = particle_t[ipart];
				send_parts[parts_page_size*num_pages_sent+part_count++] = particle_dt[ipart];
	
				for ( j = 0; j < nDim; j++ ) {
					send_parts[parts_page_size*num_pages_sent+part_count++] = particle_x[ipart][j];
				}
	
				for ( j = 0; j < nDim; j++ ) {
					send_parts[parts_page_size*num_pages_sent+part_count++] = particle_v[ipart][j];
				}

#ifdef STARFORM
				if ( particle_is_star( ipart ) ) {
					/* pack in star variables */
					send_stars[star_vars_sent++] = particle_mass[ipart];
					send_stars[star_vars_sent++] = star_initial_mass[ipart];
					send_stars[star_vars_sent++] = star_tbirth[ipart];

#ifdef ENRICH
					send_stars[star_vars_sent++] = star_metallicity_II[ipart];
#ifdef ENRICH_SNIa
					send_stars[star_vars_sent++] = star_metallicity_Ia[ipart];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
				}
#endif /* STARFORM */
	
				ipart = particle_list_next[ipart];

				if ( id_count == page_size ) {
					MPI_Isend( &send_id[num_pages_sent*page_size], page_size, MPI_INT, proc, 
						proc_pages_sent, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
					MPI_Isend( &send_parts[num_pages_sent*parts_page_size], parts_page_size, 
						MPI_DOUBLE, proc, proc_pages_sent, MPI_COMM_WORLD, 
						&send_requests[num_send_requests++] );

#ifdef STARFORM
					MPI_Isend( &send_stars[star_page_start], star_vars_sent-star_page_start,
						MPI_FLOAT, proc, proc_pages_sent, MPI_COMM_WORLD,
						&send_requests[num_send_requests++] );
					star_page_start = star_vars_sent;
#endif /* STARFORM */

					proc_pages_sent++;
					num_pages_sent++;
					id_count = 0;
					part_count = 0;
				}
			}
	
			particle_list_free( particle_list_to_send[proc] );

			if ( id_count > 0 ) {
	
				MPI_Isend( &send_id[num_pages_sent*page_size], id_count, MPI_INT, proc, proc_pages_sent, 
					MPI_COMM_WORLD, &send_requests[num_send_requests++] );
				MPI_Isend( &send_parts[num_pages_sent*parts_page_size], part_count, MPI_DOUBLE, proc, 
					proc_pages_sent, MPI_COMM_WORLD, &send_requests[num_send_requests++] );

#ifdef STARFORM
				MPI_Isend( &send_stars[star_page_start], star_vars_sent-star_page_start,
					MPI_FLOAT, proc, proc_pages_sent, MPI_COMM_WORLD,
					&send_requests[num_send_requests++] );
				star_page_start = star_vars_sent;
#endif /* STARFORM */
				num_pages_sent++;
				proc_pages_sent++;
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

#ifdef STARFORM
			MPI_Wait( &recv_stars_requests[proc], &status );
			MPI_Get_count( &status, MPI_FLOAT, &total_star_vars );
			star_vars_recv = 0;
#endif /* STARFORM */

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

#ifdef STARFORM
				if ( particle_id_is_star( recv_id[proc][i] ) ) {
					cart_assert( ipart >= 0 && ipart < num_star_particles );
					cart_assert( particle_is_star(ipart) );

					/* unpack star variables */
					particle_mass[ipart] = recv_stars[proc][star_vars_recv++];
					star_initial_mass[ipart] = recv_stars[proc][star_vars_recv++];
					star_tbirth[ipart] = recv_stars[proc][star_vars_recv++];

#ifdef ENRICH
					star_metallicity_II[ipart] = recv_stars[proc][star_vars_recv++];
#ifdef ENRICH_SNIa
					star_metallicity_Ia[ipart] = recv_stars[proc][star_vars_recv++];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
				} else {
					particle_mass[ipart] = particle_species_mass[ particle_specie( recv_id[proc][i] ) ];
				}
#else
				particle_mass[ipart] = particle_species_mass[ particle_specie( recv_id[proc][i] ) ];
#endif /* STARFORM */

				icell = cell_find_position( particle_x[ipart] );
				cart_assert( icell != -1 && cell_is_local(icell) );
				insert_particle( icell, ipart );
			}

#ifdef STARFORM
			/* ensure we unpacked all stars */
			cart_assert( star_vars_recv == total_star_vars );
#endif /* STARFORM */

			num_parts_to_recv[proc] -= id_count;

			/* if we received a full page, set up to receive a new one */
			if ( num_parts_to_recv[proc] > 0 ) {
				page_count[proc]++;
				MPI_Irecv( recv_id[proc], page_size, MPI_INT, proc, page_count[proc], 
						MPI_COMM_WORLD, &recv_id_requests[proc] );
				MPI_Irecv( recv_parts[proc], parts_page_size, MPI_DOUBLE, proc, 
						page_count[proc], MPI_COMM_WORLD, &recv_parts_requests[proc] );

#ifdef STARFORM
				MPI_Irecv( recv_stars[proc], star_page_size, MPI_FLOAT, proc, 
						page_count[proc], MPI_COMM_WORLD, &recv_stars_requests[proc] );
#endif /* STARFORM */
			} else {
				cart_free( recv_id[proc] );
				cart_free( recv_parts[proc] );

#ifdef STARFORM
				cart_free( recv_stars[proc] );
#endif /* STARFORM */
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

#ifdef STARFORM
	cart_free( send_stars );
#endif /* STARFORM */

	end_time( TRADE_PARTICLE_TIMER );
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

#ifdef STARFORM
	if ( particle_id_is_star(id) ) {
		if ( free_star_particle_list == NULL_PARTICLE ) {
			if ( next_free_star_particle >= num_star_particles ) {
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
					cart_error("Ran out of star particles, increase num_star_particles!");
				}
			} else {
				ipart = next_free_star_particle;
				next_free_star_particle++;
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
					if ( next_free_star_particle >= num_star_particles ) {
						cart_error("Ran out of particles, increase num_particles!");
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
			cart_error("Ran out of particles, increase num_particles!");
		}

		ipart = next_free_particle;
		next_free_particle++;
	} else {
		ipart = free_particle_list;
		free_particle_list = particle_list_next[free_particle_list];
	}
#endif /* STARFORM */

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

#ifdef STARFORM
	/* move star variables here */
	if ( particle_is_star(ipart_old) ) {
		cart_assert( ipart_new < num_star_particles );

		star_tbirth[ipart_new] = star_tbirth[ipart_old];
		star_initial_mass[ipart_new] = star_initial_mass[ipart_old];
		
#ifdef ENRICH
		star_metallicity_II[ipart_new] = star_metallicity_II[ipart_old];
#ifdef ENRICH_SNIa
		star_metallicity_Ia[ipart_new] = star_metallicity_Ia[ipart_old];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
	
	}
#endif /* STARFORM */

	particle_id[ipart_old] = NULL_PARTICLE;
	particle_level[ipart_old] = FREE_PARTICLE_LEVEL;
}

void particle_free( int ipart ) {
	cart_assert( ipart >= 0 && ipart < num_particles );
	cart_assert( particle_level[ipart] != FREE_PARTICLE_LEVEL );

#ifdef STARFORM
	if ( ipart < num_star_particles ) {
		particle_list_next[ipart] = free_star_particle_list;
		free_star_particle_list = ipart;

		if ( particle_is_star(ipart) ) {	
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
#endif /* STARFORM */

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

#ifdef STARFORM
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
#endif /* STARFORM */

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

int particle_specie( int id ) {
	int specie = 0;

	for ( specie = 0; specie < num_particle_species; specie++ ) {
		if ( id < particle_species_indices[specie+1] ) {
			break;
		}
	}

	return specie;
}

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)

void build_mesh() {
	int i,j;
	int level, cell;
	int total_cells_per_level[max_level-min_level+1];
	float refmin[nDim];
	float refmax[nDim];

	/* Doug (11/29/2009): necessary to properly set particle timestep 
	 * (not certain abox and auni are required here) */
	/* NG: can NOT be used here, hydro vars may not have been read yet! */
	//set_timestepping_scheme();

#ifdef REFINEMENT
	if ( spatially_limited_refinement ) {
		for ( i = 0; i < nDim; i++ ) {
			refmin[i] = num_grid+1.0;
			refmax[i] = -1.0;

			for ( j = 0; j < num_particles; j++ ) {
				if ( particle_level[j] != FREE_PARTICLE_LEVEL &&
						particle_id[j] < particle_species_indices[1] ) {
					if ( particle_x[j][i] < refmin[i] ) {
						refmin[i] = particle_x[j][i];
					}

					if ( particle_x[j][i] > refmax[i] ) {
						refmax[i] = particle_x[j][i];
					}
				}
			}

			refmin[i] = floor(refmin[i]);
			refmax[i] = ceil(refmax[i]);
		}

		MPI_Allreduce( refmin, refinement_volume_min, nDim, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
		MPI_Allreduce( refmax, refinement_volume_max, nDim, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );

		for ( i = 0; i < nDim; i++ ) {
			cart_debug("refinement_volume[%u] = %e %e", i, refinement_volume_min[i], refinement_volume_max[i] );
		}
	}
#endif /* REFINEMENT */

	build_cell_buffer();
	repair_neighbors();

#ifdef REFINEMENT
	/* do initial refinement */
	level = min_level;
	total_cells_per_level[min_level] = num_root_cells;
	while ( level < max_level && total_cells_per_level[level] > 0 ) {
		cart_debug("assigning density to level %u", level );
		assign_density(level);
		cart_debug("refining level %u, num_cells_per_level = %d", level, num_cells_per_level[level] );
		modify( level, 0 );
		cart_debug("done refining level %u, created %u new cells", 
				level, num_cells_per_level[level+1] );
		MPI_Allreduce( &num_cells_per_level[level+1], &total_cells_per_level[level+1], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
		level++;
	
		if ( local_proc_id == MASTER_NODE ) {
			cart_debug("level %u: %u cells", level, total_cells_per_level[level] );
		}

		if(total_cells_per_level[level] > 0) {
			for(j=0; j<num_particles; j++) { 
				if ( particle_level[j] != FREE_PARTICLE_LEVEL) {
					cell = cell_find_position_above_level(level,particle_x[j]);
					cart_assert(cell > -1);
					particle_level[j] = cell_level(cell);
				}
			}
		   
			load_balance();
		}
	}
	}
#else
	load_balance();
#endif /* REFINEMENT */
}

#endif /* GRAVITY || RADIATIVE_TRANSFER */

#endif /* PARTICLES */
