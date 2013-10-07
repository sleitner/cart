#include "config.h"
#if defined(STAR_FORMATION) && defined(AGN)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "agn.h"
#include "agn_step.h"
#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "halos.h"
#include "halo_finder.h"
#include "hydro.h"
#include "io.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "particle_buffer.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "sfc.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "step.h"

extern int agn_accretion_recipe;
extern int agn_feedback_recipe;
extern int agn_feedback_storage;
extern int agn_dv_on_off;

extern double eddington_factor;
extern double radiative_efficiency;
extern double feedback_efficiency;
extern double bondi_normalization;
extern double bondi_pivot_gas_density;
extern double bondi_exponent;
extern double minimum_agn_feedback_temperature;
extern double sink_particle_delta;
extern double agn_merge_radius;
extern double agn_merge_velocity;
extern double agn_seed_mass;
extern double agn_halo_mass_threshold;

double Meddfact;
double Mbondifact;
double Efbfact;
double Eagncritfact; 
double crit_bondi_gas_density;

typedef struct {
	int agn_index;
	int agn_cell;
	float sink_volume;

	int sink_cell_proc[SINK_SAMPLES_3D];
	int sink_cell_index[SINK_SAMPLES_3D];
	int sink_cell_level[SINK_SAMPLES_3D];

	float sink_density[SINK_SAMPLES_3D];
	float sink_energy[SINK_SAMPLES_3D];
	float sink_momentum[SINK_SAMPLES_3D][nDim];
} agn_sink_data;

void agn_collect_sink_values( int num_level_agn, agn_sink_data *agn_list , int level );
void agn_update_sink_values( int num_level_agn, agn_sink_data *agn_list , int level );
void agn_compute_agn_physics( int num_level_agn, agn_sink_data *agn_list , int level );
void construct_agn_root_cell_lists( int *num_agn_ret, int **num_root_agn_ret, int **agn_index_ret, int **agn_particle_list_ret );
int agn_merge_or_not( int ipart, int ipart2 );
void agn_merge_list( int ipart, int *merge_list, int merge_count );

extern double sink_sample_weights[SINK_SAMPLES_3D];
extern double sink_sample_offset[SINK_SAMPLES_3D][nDim];
extern int num_sink_samples;

void agn_feedback( int level ) {
	/* primary entry point for agn module */
	int i;
	int ipart;
	int num_level_cells, *level_cells;
	int num_level_agn;
	agn_sink_data *agn_list;
		
	double t_next = tl[level] + dtl[level]; 

	/* compute redshift-dependent factors */
	crit_bondi_gas_density = bondi_pivot_gas_density/constants->XH/units->number_density;
	Meddfact = eddington_factor*4.0*M_PI*constants->G*
		constants->mp/radiative_efficiency/constants->sigmaT/constants->c*units->time; 
	Mbondifact = 4.0*M_PI*0.52*pow(constants->G/units->length/units->length/units->length*units->mass*units->time*units->time,2.0);
	Efbfact = radiative_efficiency*feedback_efficiency*pow(constants->c/units->velocity,2.0)*cell_volume_inverse[level];

	switch ( agn_feedback_storage ) { /* agn_feedback_storage = 0 to turn off, 1 to turn on */
		case 0:
			Eagncritfact = 0.0; 
			break;
		case 1: 
			Eagncritfact = minimum_agn_feedback_temperature/units->temperature/constants->wmu/(constants->gamma-1)*cell_volume_inverse[level]; 
			/* Note: mH AND k (kB) is included in units->temperature 
             * Also, to match up in comparison to Efb/cell_mass, need to mult. by cell_volume_inverse[level] (??) */
			break;
		default: 
			cart_error("Invalid agn_feedback_storage value: %d", agn_feedback_storage );
	}

	/* select black hole particles */
	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

	num_level_agn = 0;
	for ( i = 0; i < num_level_cells; i++ ) {
		ipart = cell_particle_list[ level_cells[i] ];

		while ( ipart != NULL_PARTICLE ) {
            if ( particle_is_star(ipart) && star_particle_type[ipart] == STAR_TYPE_AGN &&
					particle_t[ipart] < t_next - 0.5*dtl[max_level] ) {	
				num_level_agn++;
			}
			ipart = particle_list_next[ipart];
		}
	}

	/* allocate agn_sink_data */
	agn_list = cart_alloc( agn_sink_data, num_level_agn );

	num_level_agn = 0;
	for ( i = 0; i < num_level_cells; i++ ) {
		ipart = cell_particle_list[ level_cells[i] ];

		while ( ipart != NULL_PARTICLE ) {
			if ( particle_is_star(ipart) && star_particle_type[ipart] == STAR_TYPE_AGN &&
					particle_t[ipart] < t_next - 0.5*dtl[max_level] ) { 
				agn_list[num_level_agn].agn_index = ipart;
				agn_list[num_level_agn++].agn_cell = cell_find_position( particle_x[ipart] );
			}
			ipart = particle_list_next[ipart];
		}
	}
	cart_free( level_cells );

	/* gather step */
	agn_collect_sink_values( num_level_agn, agn_list, level );

	/* calculation step */
	agn_compute_agn_physics( num_level_agn, agn_list, level );	

	/* update step */
	agn_update_sink_values( num_level_agn, agn_list, level );

	/* free agn_sink_data */
	cart_free( agn_list );
}

void agn_collect_sink_values( int num_level_agn, agn_sink_data *agn_list, int level ) {
	int i, j, k;
	int ipart;
	int icell;
	int proc;
	int agn_index;
	int sample;
	int sfc;
	int num_agn_cell_vars = 5;
	double r_K = 0.0; 
	double sample_x[nDim];
	int num_send_requests = 0;	
	int sync_flag[MAX_PROCS];
	int recv_sample_count[MAX_PROCS];
	int send_sample_count[MAX_PROCS];
	MPI_Request send_requests[5*MAX_PROCS];
	MPI_Request recv_sample_requests[MAX_PROCS];
	MPI_Request recv_cell_values_requests[MAX_PROCS];
	MPI_Request recv_cell_indices_requests[MAX_PROCS];
	MPI_Request recv_cell_levels_requests[MAX_PROCS];
	int *sample_agn_index[MAX_PROCS];
	int *sample_agn_point[MAX_PROCS];
	double *send_positions[MAX_PROCS];
	double *recv_positions;
	int *send_cell_indices[MAX_PROCS];
	int *recv_cell_indices[MAX_PROCS];
	int *send_cell_levels[MAX_PROCS];
	int *recv_cell_levels[MAX_PROCS];
	float *send_cell_values[MAX_PROCS];
	float *recv_cell_values[MAX_PROCS];

	/* calculate processors to communicate with */
	for ( i = 0; i < num_procs; i++ ) {
		if ( num_local_buffers[level][i] > 0 || num_remote_buffers[level][i] > 0 ) {
			sync_flag[i] = 1;
		} else {
			sync_flag[i] = 0;
		}
	}	

	for ( i = 0; i < num_procs; i++ ) {
		send_sample_count[i] = 0;
		recv_sample_count[i] = 0;

		if ( sync_flag[i] ) {
			MPI_Irecv( &recv_sample_count[i], 1, MPI_INT, i, 0, 
					mpi.comm.run, &recv_sample_requests[i] );
		} else {
			recv_sample_requests[i] = MPI_REQUEST_NULL;
		}

	}

	/* loop over black holes */
	for ( i = 0; i < num_level_agn; i++ ) {
		ipart = agn_list[i].agn_index;

		if ( SINK_RADIUS_SAMPLES == 0 ) {
			agn_list[i].sink_volume = cell_volume[level];
		} else {
			r_K = compute_r_K( ipart, agn_list[i].agn_cell );
			agn_list[i].sink_volume = pow( (SINK_ACCRETION_RADIUS+1)*r_K, 3.0 ) / SINK_SAMPLES_3D;
		}

		for ( j = 0; j < num_sink_samples; j++ ) {  
			for ( k = 0; k < nDim; k++ ) {
				sample_x[k] = ( particle_x[ipart][k] + sink_sample_offset[j][k]*r_K );
				if ( sample_x[k] >= (double)num_grid ) {
					sample_x[k] -= (double)num_grid;
				} else if ( sample_x[k] < 0.0 ) {
					sample_x[k] += (double)num_grid;
				} 
			} 

			sfc = sfc_index_position( sample_x );
			proc = processor_owner( sfc );

			if ( proc == local_proc_id ) {
				icell = cell_find_positioni_sfc( sfc, sample_x );

				agn_list[i].sink_cell_proc[j] = local_proc_id;
				agn_list[i].sink_cell_index[j] = icell;
				agn_list[i].sink_cell_level[j] = cell_level(icell);
				agn_list[i].sink_density[j] = cell_gas_density(icell);
				agn_list[i].sink_energy[j] = cell_gas_sound_speed(icell); 

				for ( k = 0; k < nDim; k++ ) {
					agn_list[i].sink_momentum[j][k] = cell_momentum(icell, k);
				}
			} else if ( proc == -1 ) {
				cart_error("Unable to find processor owner of sink position for agn %ld\n", particle_id[ipart] );
			} else {
				send_sample_count[proc]++;
				agn_list[i].sink_cell_proc[j] = proc;
			}
		}
	}

	for ( i = 0; i < num_procs; i++ ) {
		if ( sync_flag[i] ) {
			MPI_Isend( &send_sample_count[i], 1, MPI_INT, i, 0,
					mpi.comm.run, &send_requests[num_send_requests++] );
			if ( send_sample_count[i] > 0 ) {
				sample_agn_index[i] = cart_alloc( int, send_sample_count[i] );
				sample_agn_point[i] = cart_alloc( int, send_sample_count[i] );
				send_positions[i] = cart_alloc( double, nDim*send_sample_count[i] );
			}
		} else {
			cart_assert( send_sample_count[i] == 0 );
		}

		send_sample_count[i] = 0;
	}

	for ( i = 0; i < num_level_agn; i++ ) {  
		ipart = agn_list[i].agn_index;

		if ( SINK_RADIUS_SAMPLES == 0 ) {
			r_K = 0.0;
		} else {
			r_K = compute_r_K( ipart, agn_list[i].agn_cell );
		}

		for ( j = 0; j < num_sink_samples; j++ ) {
			if ( agn_list[i].sink_cell_proc[j] != local_proc_id ) {
				for ( k = 0; k < nDim; k++ ) {
					sample_x[k] = ( particle_x[ipart][k] + sink_sample_offset[j][k]*r_K );
					if ( sample_x[k] >= (double)num_grid ) {
						sample_x[k] -= (double)num_grid;
					} else if ( sample_x[k] < 0.0 ) {
						sample_x[k] += (double)num_grid;
					}
				}

				proc = agn_list[i].sink_cell_proc[j];
				cart_assert( sync_flag[proc] );

				sample_agn_index[proc][send_sample_count[proc]] = i;
				sample_agn_point[proc][send_sample_count[proc]] = j;

				for ( k = 0; k < nDim; k++ ) {
					send_positions[proc][nDim*send_sample_count[proc]+k] = sample_x[k];
				}

				send_sample_count[proc]++;
			}
		}
	}

	for ( i = 0; i < num_procs; i++ ) {
		if ( send_sample_count[i] > 0 ) {
			MPI_Isend( send_positions[i], nDim*send_sample_count[i], MPI_DOUBLE,
					i, 0, mpi.comm.run, &send_requests[num_send_requests++] );
			recv_cell_indices[i] = cart_alloc( int, send_sample_count[i] );
			recv_cell_levels[i] = cart_alloc( int, send_sample_count[i] );
			recv_cell_values[i] = cart_alloc(float, num_agn_cell_vars*send_sample_count[i] );

			MPI_Irecv( recv_cell_indices[i], send_sample_count[i], MPI_INT,
					i, 1, mpi.comm.run, &recv_cell_indices_requests[i] );
			MPI_Irecv( recv_cell_levels[i], send_sample_count[i], MPI_INT,
					i, 2, mpi.comm.run, &recv_cell_levels_requests[i] );
			MPI_Irecv( recv_cell_values[i], num_agn_cell_vars*send_sample_count[i], 
					MPI_FLOAT, i, 1, mpi.comm.run, 
					&recv_cell_values_requests[i] );
		} else {
			recv_cell_indices_requests[i] = MPI_REQUEST_NULL;
			recv_cell_levels_requests[i] = MPI_REQUEST_NULL;
			recv_cell_values_requests[i] = MPI_REQUEST_NULL;
		}
	}

	do {
		MPI_Waitany( num_procs, recv_sample_requests, &i, MPI_STATUS_IGNORE );

		if ( i != MPI_UNDEFINED ) {
			if ( recv_sample_count[i] > 0 ) {
				recv_positions = cart_alloc( double, nDim*recv_sample_count[i] );

				MPI_Recv( recv_positions, nDim*recv_sample_count[i], MPI_DOUBLE,
						i, 0, mpi.comm.run, MPI_STATUS_IGNORE );

				send_cell_indices[i] = cart_alloc( int, recv_sample_count[i] );
				send_cell_levels[i] = cart_alloc( int, recv_sample_count[i] );
				send_cell_values[i] = cart_alloc( float, num_agn_cell_vars*recv_sample_count[i] );

				/* process received sample requests */
				for ( j = 0; j < recv_sample_count[i]; j++ ) {
					icell = cell_find_position( &recv_positions[nDim*j] ); 
					cart_assert( icell != NULL_OCT );

					send_cell_indices[i][j] = icell;
					send_cell_levels[i][j] = cell_level(icell);
					send_cell_values[i][num_agn_cell_vars*j] = cell_gas_density(icell);
					send_cell_values[i][num_agn_cell_vars*j+1] = cell_gas_sound_speed(icell);
					send_cell_values[i][num_agn_cell_vars*j+2] = cell_momentum(icell,0);
					send_cell_values[i][num_agn_cell_vars*j+3] = cell_momentum(icell,1);
					send_cell_values[i][num_agn_cell_vars*j+4] = cell_momentum(icell,2);
				}

				MPI_Isend( send_cell_indices[i], recv_sample_count[i], MPI_INT,
						i, 1, mpi.comm.run, &send_requests[num_send_requests++] );
				MPI_Isend( send_cell_levels[i], recv_sample_count[i], MPI_INT,
						i, 2, mpi.comm.run, &send_requests[num_send_requests++] );
				MPI_Isend( send_cell_values[i], num_agn_cell_vars*recv_sample_count[i],
						MPI_FLOAT, i, 1, mpi.comm.run, &send_requests[num_send_requests++] );
				cart_free( recv_positions );
			}
		}
	} while ( i != MPI_UNDEFINED );

	do {
		MPI_Waitany( num_procs, recv_cell_indices_requests, &i, MPI_STATUS_IGNORE );

		if ( i != MPI_UNDEFINED ) {
			MPI_Wait( &recv_cell_values_requests[i], MPI_STATUS_IGNORE );
			MPI_Wait( &recv_cell_levels_requests[i], MPI_STATUS_IGNORE );

			for ( j = 0; j < send_sample_count[i]; j++ ) {
				agn_index = sample_agn_index[i][j];
				sample = sample_agn_point[i][j];

				cart_assert( sample >= 0 && sample < num_sink_samples );

				agn_list[agn_index].sink_cell_index[sample] = recv_cell_indices[i][j];
				cart_assert( recv_cell_levels[i][j] >= min_level && recv_cell_levels[i][j] <= max_level );
				agn_list[agn_index].sink_cell_level[sample] = recv_cell_levels[i][j];

				agn_list[agn_index].sink_density[sample] = recv_cell_values[i][num_agn_cell_vars*j];
				agn_list[agn_index].sink_energy[sample] = recv_cell_values[i][num_agn_cell_vars*j+1];
				agn_list[agn_index].sink_momentum[sample][0] = recv_cell_values[i][num_agn_cell_vars*j+2];
				agn_list[agn_index].sink_momentum[sample][1] = recv_cell_values[i][num_agn_cell_vars*j+3];
				agn_list[agn_index].sink_momentum[sample][2] = recv_cell_values[i][num_agn_cell_vars*j+4];
			}
			cart_free( recv_cell_indices[i] );
			cart_free( recv_cell_levels[i] );
			cart_free( recv_cell_values[i] );
			cart_free( sample_agn_index[i] );
			cart_free( sample_agn_point[i] );
		}
	} while ( i != MPI_UNDEFINED );

	/* wait for all sends to complete before freeing buffers */
	MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	
	for ( i = 0; i < num_procs; i++ ) {
		if ( send_sample_count[i] > 0 ) {
			cart_free( send_positions[i] );
		}

		if ( recv_sample_count[i] > 0 ) {
			cart_free( send_cell_indices[i] );
			cart_free( send_cell_levels[i] );
			cart_free( send_cell_values[i] );
		}
	}
}

void agn_update_sink_values( int num_level_agn, agn_sink_data *agn_list, int level ) {
	int i, j, k;
	int ipart, icell;
	int proc, sample;
	int num_agn_cell_vars = 5;
	int num_send_requests = 0;	
	double new_density, density_fraction;
	int sync_flag[MAX_PROCS];
	int recv_sample_count[MAX_PROCS];
	int send_sample_count[MAX_PROCS];
	MPI_Request send_requests[3*MAX_PROCS];
	MPI_Request recv_sample_requests[MAX_PROCS];
	int *send_cell_indices[MAX_PROCS];
	float *send_cell_deltas[MAX_PROCS];
	int *recv_cell_indices[MAX_PROCS];
	float *recv_cell_deltas[MAX_PROCS];

	/* calculate processors to communicate with */
	for ( i = 0; i < num_procs; i++ ) {
		if ( num_local_buffers[level][i] > 0 || num_remote_buffers[level][i] > 0 ) {
			sync_flag[i] = 1;
		} else {
			sync_flag[i] = 0;
		}
	}

	/* set up receives */
	for ( i = 0; i < num_procs; i++ ) { 
		send_sample_count[i] = 0;
		recv_sample_count[i] = 0;

		if ( sync_flag[i] ) {
			MPI_Irecv( &recv_sample_count[i], 1, MPI_INT, i, 0, 
					mpi.comm.run, &recv_sample_requests[i] );
		} else {
			recv_sample_requests[i] = MPI_REQUEST_NULL;
		}
	}

	for ( i = 0; i < num_level_agn; i++ ) {
        for ( j = 0; j < num_sink_samples; j++ ) {
            if ( agn_list[i].sink_cell_proc[j] != local_proc_id ) {
                send_sample_count[ agn_list[i].sink_cell_proc[j] ]++;
            }
        }
    }  

	for ( i = 0; i < num_procs; i++ ) {
		if ( sync_flag[i] ) {
			MPI_Isend( &send_sample_count[i], 1, MPI_INT, i, 0,
					mpi.comm.run, &send_requests[num_send_requests++] );

			if ( send_sample_count[i] > 0 ) {
				send_cell_indices[i] = cart_alloc( int, send_sample_count[i] );
				send_cell_deltas[i] = cart_alloc( float, num_agn_cell_vars*send_sample_count[i] );
			} 
		} else {
			cart_assert( send_sample_count[i] == 0 );
		}

		send_sample_count[i] = 0;
	}

	for ( i = 0; i < num_level_agn; i++ ) {
		ipart = agn_list[i].agn_index;

		for ( j = 0; j < num_sink_samples; j++ ) {
			if ( agn_list[i].sink_cell_proc[j] != local_proc_id ) {
				proc = agn_list[i].sink_cell_proc[j];
				sample = send_sample_count[proc];
				send_cell_indices[proc][sample] = agn_list[i].sink_cell_index[j]; 
				send_cell_deltas[proc][num_agn_cell_vars*sample] = agn_list[i].sink_density[j];
				send_cell_deltas[proc][num_agn_cell_vars*sample+1] = agn_list[i].sink_energy[j];
				send_cell_deltas[proc][num_agn_cell_vars*sample+2] = agn_list[i].sink_momentum[j][0];
				send_cell_deltas[proc][num_agn_cell_vars*sample+3] = agn_list[i].sink_momentum[j][1];
				send_cell_deltas[proc][num_agn_cell_vars*sample+4] = agn_list[i].sink_momentum[j][2];
				send_sample_count[proc]++;
			} else {
				icell = agn_list[i].sink_cell_index[j];
				new_density = cell_gas_density(icell) - 
					agn_list[i].sink_density[j] * cell_volume_inverse[ cell_level(icell) ];
				density_fraction = new_density / cell_gas_density(icell);

				cell_gas_density(icell) = new_density;
				cell_gas_energy(icell) *= density_fraction;
				cell_gas_internal_energy(icell) *= density_fraction;

				for ( k = 0; k < nDim; k++ ) {
					cell_momentum(icell,k) = density_fraction*cell_momentum(icell,k) + agn_list[i].sink_momentum[j][k];
				}

				for ( k = 0; k < num_chem_species; k++ ) {                                                                                                                                  
					cell_advected_variable(icell,k) *= density_fraction;
				}

				/* energy due to feedback */
				cell_gas_energy(icell) += agn_list[i].sink_energy[j];
				cell_gas_internal_energy(icell) += agn_list[i].sink_energy[j];
				cell_gas_pressure(icell) = MAX( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );
			}
		}
	}

	for ( i = 0; i < num_procs; i++ ) {
		if ( send_sample_count[i] > 0 ) {
			MPI_Isend( send_cell_indices[i], send_sample_count[i], MPI_INT, i, 
					1, mpi.comm.run, &send_requests[num_send_requests++] );
			MPI_Isend( send_cell_deltas[i], num_agn_cell_vars*send_sample_count[i], 
					MPI_FLOAT, i, 1, mpi.comm.run, &send_requests[num_send_requests++] );
		}
	}

	do {
		MPI_Waitany( num_procs, recv_sample_requests, &i, MPI_STATUS_IGNORE );
		if ( i != MPI_UNDEFINED ) {
			if ( recv_sample_count[i] > 0 ) {
				recv_cell_indices[i] = cart_alloc( int, recv_sample_count[i] );
				recv_cell_deltas[i] = cart_alloc( float, num_agn_cell_vars*recv_sample_count[i] );

				MPI_Recv( recv_cell_indices[i], recv_sample_count[i], MPI_INT, i, 
						1, mpi.comm.run, MPI_STATUS_IGNORE );
				MPI_Recv( recv_cell_deltas[i], num_agn_cell_vars*recv_sample_count[i], 
						MPI_FLOAT, i, 1, mpi.comm.run, MPI_STATUS_IGNORE );

				/* process received sample requests */
				for ( j = 0; j < recv_sample_count[i]; j++ ) {
					icell = recv_cell_indices[i][j];
					new_density = cell_gas_density(icell) - 
						recv_cell_deltas[i][j*num_agn_cell_vars] * cell_volume_inverse[ cell_level(icell) ];
					/* Note that sink_density here has units of mass (replaced list quantities) */
					density_fraction = new_density / cell_gas_density(icell);

					cell_gas_density(icell) = new_density;
					cell_gas_energy(icell) *= density_fraction;
					cell_gas_internal_energy(icell) *= density_fraction;

					for ( k = 0; k < nDim; k++ ) {
						cell_momentum(icell,k) = density_fraction*cell_momentum(icell,k) + recv_cell_deltas[i][j*num_agn_cell_vars+2+k];
					}

					for ( k = 0; k < num_chem_species; k++ ) {                                                                                                                                  
						cell_advected_variable(icell,k) *= density_fraction;
					}

					/* energy due to feedback */
					cell_gas_energy(icell) += recv_cell_deltas[i][j*num_agn_cell_vars+1];
					cell_gas_internal_energy(icell) += recv_cell_deltas[i][j*num_agn_cell_vars+1];
					cell_gas_pressure(icell) = MAX( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );
				}

				cart_free( recv_cell_indices[i] ); 
				cart_free( recv_cell_deltas[i] );
			}
		}
	} while ( i != MPI_UNDEFINED );

	/* wait for all sends to complete before freeing buffers */
	MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	for ( i = 0; i < num_procs; i++ ) {
		if ( send_sample_count[i] > 0 ) {
			cart_free( send_cell_indices[i] );
			cart_free( send_cell_deltas[i] );
		}
	}
}


void agn_compute_agn_physics( int num_level_agn, agn_sink_data *agn_list, int level  ) {
	int i,j,k;
	int ipart;
	double sink_averaged_cs, sink_averaged_density; 
	double sink_averaged_gas_velocity[nDim];
	double change_in_mass[SINK_SAMPLES_3D], change_in_momentum[SINK_SAMPLES_3D][nDim]; 
	double Macc, newMacc, Medd, new_momentum[nDim], new_particle_mass, Mbondi, Mbondibetafact, Mbh, Efb, virtual_sample_volume, sink_mass;
	double dv = 0.0; 
	double delta_t;

	double t_next = tl[level] + dtl[level]; 

#ifdef COSMOLOGY 
#pragma omp parallel for default(none), private( i, j, k, ipart, Mbh, virtual_sample_volume, delta_t, sink_mass, sink_averaged_cs, sink_averaged_density, sink_averaged_gas_velocity, new_momentum, Medd, Macc, dv, Mbondibetafact, Mbondi, change_in_mass, newMacc, change_in_momentum, new_particle_mass, Efb ), shared( num_level_agn,  agn_list, t_next, num_sink_samples, sink_sample_weights, Meddfact, Mbondifact, Efbfact, Eagncritfact, agn_accretion_recipe, particle_x, particle_v, particle_t, units, constants, particle_mass, star_metallicity_II, cell_vars, tl, abox, level, crit_bondi_gas_density, agn_dv_on_off, particle_id, bondi_normalization, bondi_exponent, radiative_efficiency )
#else
#pragma omp parallel for default(none), private( i, j, k, ipart, Mbh, virtual_sample_volume, delta_t, sink_mass, sink_averaged_cs, sink_averaged_density, sink_averaged_gas_velocity, new_momentum, Medd, Macc, dv, Mbondibetafact, Mbondi, change_in_mass, newMacc, change_in_momentum, new_particle_mass, Efb ), shared( num_level_agn,  agn_list, t_next, num_sink_samples, sink_sample_weights, Meddfact, Mbondifact, Efbfact, Eagncritfact, agn_accretion_recipe, particle_x, particle_v, particle_t, units, constants, particle_mass, star_metallicity_II, cell_vars, tl, level, crit_bondi_gas_density, agn_dv_on_off, particle_id, bondi_normalization, bondi_exponent, radiative_efficiency )
#endif
	/* loop over black holes */
	for ( i=0; i<num_level_agn; i++ ) {
		ipart = agn_list[i].agn_index;
		cart_assert( ipart >= 0 && ipart < num_particles );

		Mbh = particle_mass[ipart];
		virtual_sample_volume = agn_list[i].sink_volume;
		delta_t = t_next - particle_t[ipart];

		sink_mass = 0.0;
		sink_averaged_cs = 0.0;
		sink_averaged_density = 0.0;
		for ( j=0; j<nDim; j++ ) {
			sink_averaged_gas_velocity[j] = 0.0;
		}

		newMacc = 0.0;  
		for ( j=0; j<nDim; j++ ) {
			new_momentum[j] = 0.0;
		}

		for ( j=0; j<num_sink_samples; j++ ) {
			/* Note: we don't need the following 6 lines if we're doing agn_accretion_recipe = 0 */
			sink_averaged_cs += agn_list[i].sink_energy[j]*sink_sample_weights[j];
			sink_averaged_density += agn_list[i].sink_density[j] * sink_sample_weights[j];
			for ( k=0; k<nDim; k++ ) {
				sink_averaged_gas_velocity[k] += sink_sample_weights[j] * ( agn_list[i].sink_momentum[j][k]/agn_list[i].sink_density[j] );
			}

			sink_mass += virtual_sample_volume * agn_list[i].sink_density[j];
		}

		Medd = Meddfact*Mbh*delta_t;

		switch ( agn_accretion_recipe ) {
			case 0: /* Pure Eddington accretion */
				Macc = Medd;
				break;
			case 1: /* Bondi accretion (constant alpha) */
				dv = (particle_v[ipart][0]-sink_averaged_gas_velocity[0])*(particle_v[ipart][0]-sink_averaged_gas_velocity[0])+
					(particle_v[ipart][1]-sink_averaged_gas_velocity[1])*(particle_v[ipart][1]-sink_averaged_gas_velocity[1])+
					(particle_v[ipart][2]-sink_averaged_gas_velocity[2])*(particle_v[ipart][2]-sink_averaged_gas_velocity[2]);

				Mbondi = bondi_normalization*Mbondifact * Mbh*Mbh * sink_averaged_density * 
							pow( agn_dv_on_off*dv + sink_averaged_cs*sink_averaged_cs, -1.5 )*delta_t; 
				Macc = MIN( Medd, Mbondi );
				break;
			case 2: /* Bondi accretion (constant beta) */
				dv = (particle_v[ipart][0]-sink_averaged_gas_velocity[0])*(particle_v[ipart][0]-sink_averaged_gas_velocity[0])+
					(particle_v[ipart][1]-sink_averaged_gas_velocity[1])*(particle_v[ipart][1]-sink_averaged_gas_velocity[1])+
					(particle_v[ipart][2]-sink_averaged_gas_velocity[2])*(particle_v[ipart][2]-sink_averaged_gas_velocity[2]);

				if ( sink_averaged_density > crit_bondi_gas_density ) {
					Mbondibetafact = pow( sink_averaged_density / crit_bondi_gas_density, bondi_exponent ) * Mbondifact; 
				} else {
					Mbondibetafact = Mbondifact;
				}
				Mbondi = Mbondibetafact * Mbh*Mbh * sink_averaged_density * 
							pow( agn_dv_on_off*dv + sink_averaged_cs*sink_averaged_cs, -1.5 )*delta_t;
				Macc = MIN( Medd, Mbondi );
				break;
			default :
				cart_error("Invalid agn_accretion_recipe in agn_compute_agn_physics");
		}

		/* Calculate change in energy (due to feedback), change_in_mass, change_in_momentum */
		for ( j=0; j<num_sink_samples; j++ ) {
			cart_assert( agn_list[i].sink_cell_level[j] >= 0 && agn_list[i].sink_cell_level[j] <= max_level );
			change_in_mass[j] = MIN( Macc * sink_sample_weights[j] , 
									.667 * agn_list[i].sink_density[j] * 
									cell_volume[agn_list[i].sink_cell_level[j]] * sink_sample_weights[j] ); 
			newMacc += change_in_mass[j];

			for ( k=0; k<nDim; k++ ) {
				change_in_momentum[j][k] = ( change_in_mass[j] / agn_list[i].sink_density[j] ) * agn_list[i].sink_momentum[j][k]; 
				new_momentum[k] += change_in_momentum[j][k];
			}
		}

		/* Calculate new velocity and mass for black hole from accreted gas mass */
		new_particle_mass = particle_mass[ipart] + ( 1-radiative_efficiency ) * newMacc;

		for ( j=0; j<nDim; j++ ) {
			particle_v[ipart][j] = ( particle_v[ipart][j] * particle_mass[ipart] + new_momentum[j] ) / new_particle_mass;
		}

		particle_mass[ipart] = new_particle_mass;

		/* Compute feedback energy and change in energy */

		Efb = Efbfact*newMacc;
#ifndef ENRICHMENT 
#error "Enrichment required for AGN stored feedback energy"
#endif /* ENRICHMENT */
		star_metallicity_II[ipart] += Efb;  /* Add feedback energy from this timestep to stored energy */

		if ( star_metallicity_II[ipart] < Eagncritfact*sink_mass ) {  
			/* If we do not make the critical energy with feedback energy from this timestep, continue to store -  no effect on environment */
			Efb = 0.0;
		} else {  
			/* If we exceed the critical energy, dump the stored energy into feedback energy - environment will heat */
			Efb = star_metallicity_II[ipart]; 
			star_metallicity_II[ipart] = 0.0;
		}  

		/* Reuse agn list for the change values */
		for ( j=0; j<num_sink_samples; j++ ) {
			agn_list[i].sink_energy[j] = sink_sample_weights[j] * Efb;
			agn_list[i].sink_density[j] = change_in_mass[j];
			for ( k=0; k<nDim; k++ ) {
				agn_list[i].sink_momentum[j][k] = change_in_momentum[j][k];
			}
		}
	}
}

void construct_agn_root_cell_lists( int *num_agn_ret, int **num_root_agn_ret, int **agn_index_ret, int **agn_particle_list_ret ) {
	int icell, ipart;
	int num_root;
	int num_agn;
	int *num_root_agn;
	int *agn_index;
	int *agn_particle_list;

	/* allocate new lists to hold AGN */
	num_root = num_cells_per_level[min_level]+num_buffer_cells[min_level];
	num_root_agn = cart_alloc(int, num_root);
	agn_index = cart_alloc(int, num_root);

	for ( icell = 0; icell < num_root; icell++ ) {
		num_root_agn[icell] = 0;
	}

	/* collect agn into root cells */
	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL &&
				particle_is_star(ipart) && star_particle_type[ipart] == STAR_TYPE_AGN ) {
			icell = cell_find_position_level( min_level, particle_x[ipart] );
			num_root_agn[icell]++;
		}
	}

	num_agn = 0;
	for ( icell = 0; icell < num_root; icell++ ) {
		agn_index[icell] = num_agn;
		num_agn += num_root_agn[icell];
		num_root_agn[icell] = 0;
	}

	agn_particle_list = cart_alloc(int, num_agn);

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL &&
				particle_is_star(ipart) && star_particle_type[ipart] == STAR_TYPE_AGN ) {
			icell = cell_find_position_level( min_level, particle_x[ipart] );
			agn_particle_list[ agn_index[icell] + num_root_agn[icell]++ ] = ipart;
		}
	}

	*num_agn_ret = num_agn;
	*num_root_agn_ret = num_root_agn;
	*agn_index_ret = agn_index;
	*agn_particle_list_ret = agn_particle_list;
}

    
void agn_find_mergers() {
	int i, j;
	int dx, di, dj, dk;
	int coords[nDim];
	int coords2[nDim];
	int icell, ipart, ipart2;
	int num_agn;
	int agn_merger_count;
	int *agn_particle_list;
	int *agn_index;
	int *num_root_agn;
	int *order;

#define MAX_AGN_MERGERS		10
	int merge_list[MAX_AGN_MERGERS];

	/* create buffer of agn particles */
	build_particle_buffer( num_particle_species - 1, STAR_TYPE_AGN );

	/* collect root-lists of agn particles for efficient spatial search */
	construct_agn_root_cell_lists( &num_agn, &num_root_agn, &agn_index, &agn_particle_list );

	/* sort by mass */
	order = cart_alloc( int, num_agn );
    for ( i = 0; i < num_agn; i++ ) {
        order[i] = agn_particle_list[i];
    }

	qsort( order, num_agn, sizeof(int), compare_particle_mass );
	dx = (int)( agn_merge_radius * cell_size[max_level] / cell_size[min_level]) + 1;

	for ( i = num_agn-1; i >= 0; i-- ) {
		ipart = order[i];

		/* skip over AGN that have been removed due to merging with more massive AGN */
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL ) {
			cart_debug("finding mergers for agn id = %ld, Mbh = %e Msun", 
				particle_id[ipart], particle_mass[ipart] * units->mass / constants->Msun );

			agn_merger_count = 0;

			for ( j = 0; j < nDim; j++ ) {
				coords[j] = (int)(particle_x[ipart][j]);
			}

			/* not necessary since we restrict dx to 0.5 cell_size[min_level], but generalizable */
			for ( di = -dx; di <= dx; di++ ) {
				coords2[0] = (coords[0] + di) % num_grid;
				for ( dj = -dx; dj <= dx; dj++ ) {
					coords2[1] = (coords[1] + dj) % num_grid;
					for ( dk = -dx; dk <= dx; dk++ ) {
						coords2[2] = (coords[2] + dk) % num_grid;

						icell = root_cell_location( sfc_index( coords2 ) );

						if ( icell != -1 ) {
							for ( j = agn_index[icell]; j < agn_index[icell]+num_root_agn[icell]; j++ ) {
								ipart2 = agn_particle_list[j];

								/* compare the two black holes (skipping previously merged black holes) */
								if ( ipart != ipart2 &&
										particle_level[ipart2] != FREE_PARTICLE_LEVEL &&
										particle_mass[ipart2] <= particle_mass[ipart] ) {

									if ( agn_merge_or_not( ipart, ipart2 ) ) {
										cart_assert( agn_merger_count < MAX_AGN_MERGERS );
										merge_list[agn_merger_count++] = ipart2;
									}	
								}
							}
						}
					}
				}
			}

			/* merge black hole list */
			if ( agn_merger_count > 0 ) {
				agn_merge_list( ipart, merge_list, agn_merger_count );
			}
		}
	}

	cart_free( num_root_agn );
	cart_free( agn_index );
	cart_free( order );
	cart_free( agn_particle_list );

	destroy_particle_buffer();

#undef MAX_AGN_MERGERS
}

void agn_merge_list( int ipart, int *merge_list, int merge_count ) {
	int i, j;
	int old_cell, new_cell;

	old_cell = cell_find_position( particle_x[ipart] );

	for ( j = 0; j < nDim; j++ ) {
		particle_x[ipart][j] *= particle_mass[ipart];
		particle_v[ipart][j] *= particle_mass[ipart];	
	}

	/* Calculate the merged mass, center of mass, and mass weighted velocity */
	for ( i = 0; i < merge_count; i++ ) {
		particle_mass[ipart] += particle_mass[merge_list[i]];

		/* Find the center of mass and mass weighted velocity */
		for ( j = 0; j < nDim; j++ ) {
			particle_x[ipart][j] += particle_x[ merge_list[i] ][j] * particle_mass[ merge_list[i] ];
			particle_v[ipart][j] += particle_v[ merge_list[i] ][j] * particle_mass[ merge_list[i] ];
		}
	}

	/* Normalize the particle position and velocity */
	for ( j = 0; j < nDim; j++ ) {
		particle_x[ipart][j] /= particle_mass[ipart];
		particle_v[ipart][j] /= particle_mass[ipart];
	}

	/* update linked list of agn particle */
	new_cell = cell_find_position( particle_x[ipart] );
	delete_particle( old_cell, ipart );
	insert_particle( new_cell, ipart );

	/* Delete all the merged particles, keeping the largest one */
	for ( i = 0; i < merge_count; i++ ) {
		delete_particle( cell_find_position( particle_x[merge_list[i]] ), merge_list[i] );
		particle_free( merge_list[i] ) ;
	}
}

int agn_merge_or_not( int ipart, int ipart2 ) {
	int i,j;
	int nsteps = 1 << max_level; /* 2^max_level */
	double max_merging_distance;  
	double net_distance;
	double net_velocity;
	double extrap_dt = dtl[max_level];
	double x1new[nDim];
	double x2new[nDim];
	double dv1[nDim];
	double dv2[nDim];

	max_merging_distance = agn_merge_radius * cell_size[max_level];
	net_velocity = sqrt( (particle_v[ipart][0]-particle_v[ipart2][0]) *
			(particle_v[ipart][0]-particle_v[ipart2][0]) + 
			(particle_v[ipart][1]-particle_v[ipart2][1]) *
			(particle_v[ipart][1]-particle_v[ipart2][1]) + 
			(particle_v[ipart][2]-particle_v[ipart2][2]) *
			(particle_v[ipart][2]-particle_v[ipart2][2]) );

	if (net_velocity*units->velocity/constants->kms < agn_merge_velocity ) {
		for ( j = 0; j < nDim; j++ ) {
			x1new[j] = particle_x[ipart][j];
			x2new[j] = particle_x[ipart2][j];
			dv1[j] = particle_v[ipart][j] * extrap_dt; 
			dv2[j] = particle_v[ipart2][j] * extrap_dt;
		}

		for ( i = 0; i < nsteps; i++ ) {
			for ( j = 0; j < nDim; j++ ) {
				x1new[j] += dv1[j];
				x2new[j] += dv2[j];

				if ( x1new[j] >= num_grid ) {
					x1new[j] -= num_grid;
				} else if ( x1new[j] <= num_grid ) {
					x1new[j] += num_grid;
				}
				if ( x2new[j] >= num_grid ) {
					x2new[j] -= num_grid;
				} else if ( x2new[j] <= num_grid ) {
					x2new[j] += num_grid;
				}       
			}

			net_distance = compute_distance_periodic( x1new, x2new );
			if ( net_distance < max_merging_distance ) {
				return 1;
			}     
		}
	}

	return 0;
}

void agn_seed( halo_list *list ) {
	int i, j, k, m;
	int ih;
	int icell, ipart;
	int num_agn;
	int *num_root_agn;
	int *agn_index;
	int *agn_particle_list;
	int coords[nDim];
	int level;
	double r;
	int seed_flag, global_seed_flag;
	halo_list *halos;
	halo *h;
	double halo_mass_threshold;
	double seed_mass;

	/* find halos */
	if ( list == NULL ) {
		halos = find_halos();
	} else {
		halos = list;
	}

	/* collect root-lists of agn particles for efficient spatial search */
	construct_agn_root_cell_lists( &num_agn, &num_root_agn, &agn_index, &agn_particle_list );

	halo_mass_threshold = agn_halo_mass_threshold*constants->Msun/units->mass;
	seed_mass = agn_seed_mass*constants->Msun/units->mass;

	for ( ih = 0; ih < halos->num_halos; ih++ ) {
		h = &halos->list[ih];

		if ( h->mvir >= halo_mass_threshold ) {
			seed_flag = 1;

			/* check for existing AGN within rvir */
			for ( i = (int)floor(h->pos[0]-h->rvir); i <= (int)(h->pos[0]+h->rvir); i++ ) {
				coords[0] = ( i + num_grid ) % num_grid;
				for ( j = (int)floor(h->pos[1]-h->rvir); j <= (int)(h->pos[1]+h->rvir); j++ ) {
					coords[1] = ( j + num_grid ) % num_grid;
					for ( k = (int)floor(h->pos[2]-h->rvir); k <= (int)(h->pos[2]+h->rvir); k++ ) {
						coords[2] = ( k + num_grid ) % num_grid;
						icell = root_cell_location( sfc_index( coords ) );
						if ( icell != NULL_OCT ) {
							for ( m = agn_index[icell]; m < agn_index[icell]+num_root_agn[icell]; m++ ) {
								ipart = agn_particle_list[m];
								r = compute_distance_periodic( particle_x[ipart], h->pos );
								if ( r < h->rvir ) {
									/* no need to seed, halo already contains a black hole */
									seed_flag = 0;
									break;
								}
							}
						}
					}
				}
			}

			MPI_Allreduce( &seed_flag, &global_seed_flag, 1, MPI_INT, MPI_MIN, mpi.comm.run );

			if ( global_seed_flag ) {
				icell = cell_find_position( h->pos );

				if ( icell != NULL_OCT && cell_is_local(icell) ) {
					ipart = create_star_particle( icell, seed_mass, dtl[cell_level(icell)], STAR_TYPE_AGN );
					cart_debug("black hole with mass %e seeded into halo of mass %e", 
							particle_mass[ipart]*units->mass/constants->Msun, 
							h->mvir*units->mass/constants->Msun );
				}
			}
		}
	}

	/* update buffered hydrodynamic variables */
	for ( level = max_level; level >= min_level; level-- ) {
		if ( level < max_level ) {
			hydro_split_update(level);
		}
		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	}   

	cart_free( num_root_agn );
	cart_free( agn_index );
	cart_free( agn_particle_list );

	/* destroy halo list */
	if ( list == NULL ) {
		destroy_halo_list(halos);
	}
}


#endif /* STAR_FORMATION && AGN */
