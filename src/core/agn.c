#include "config.h"
#if defined(STARFORM) && defined(AGN)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "io.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "particle_buffer.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "tree.h"
#include "units.h"
#include "agn.h"

/* agn model parameters */
int agn_accretion_recipe = 1;       /* Set to 0 for pure eddington accretion, 1 for Bondi accretion (constant-alpha), 2 for Bondi accretion (constant-beta) */
int agn_feedback_recipe = 0;
int agn_feedback_storage = 1;       /* Set to 0 for no storage, 1 for storage via Booth and Schaye (2009) Eqn (7) */
int agn_dv_on_off = 1;              /* Set dv to 0 in denominator of accretion rate or not (for Bondi-type accretion) */
int agn_merge_radius_type = 0;      /* 0 for multiples of resolution size, 1 for comoving units */

struct {
	double eddington_factor;        /* normalization of Eddington accretion (1 for standard value) */
	double bondi_normalization;		/* alpha normalization of standard Bondi implementations */
	double radiative_efficiency;    /* eta or eps_r */
	double feedback_efficiency;	    /* eps_f */
	double critical_SF_gas_density; /* n_H^* */  
	double bondi_exponent;          /* Beta for constant-Beta (density dependent) model of Bondi implementations */
	double minimum_feedback_temperature;   /* Tmin for critical energy in feedback */  
	double sink_particle_delta;     /* This determines maximum r_acc=sink_particle_delta*resolution_size */
	double agn_merge_radius;        /* maximum distance to merge two black hole particles */
	double agn_merge_velocity;      /* maximum velocity to merge in km/s */
}

/* default parameter values */
agn_parameters = { 1.0, 100.0, 0.1, 0.05 , 0.1 , 2 , 1e8, 4.0, 2, 1000 }; 

void config_init_agn() {

	control_parameter_add(control_parameter_double,&agn_parameters.eddington_factor,
			"agn:eddington_factor",
			"Normalization of Eddington accretion rate (typically a small number 1-few)");
	control_parameter_add2(control_parameter_double,&agn_parameters.bondi_normalization,
			"agn:bondi_normalization","agn:alpha",
			"Normalization of Bondi-Hoyle accretion rate, used to compensate for not resolving r_bondi" );
	control_parameter_add3(control_parameter_double,&agn_parameters.radiative_efficiency,
			"agn:radiative_efficiency","agn:eps_r","agn:eta",
			"Radiative efficiency, L = \\eta \\dot{M} c^2" );
	control_parameter_add2(control_parameter_double,&agn_parameters.feedback_efficiency,
			"agn:feedback_efficiency","agn:eps_f",
			"Feedback efficiency, \\dot{E} = eps_f L = \\eps_f \\eps_r \\dot{M} c^2" );
	control_parameter_add(control_parameter_int, &agn_accretion_recipe,
			      "agn:accretion_recipe", "Accretion recipe for AGN: 0 = pure Eddington, 1 = Bondi Constant alpha model in single cell, 2= Bondi Constant beta model in single cell" ); 
	control_parameter_add2(control_parameter_double,&agn_parameters.critical_SF_gas_density,"agn:critical_SF_gas_density","agn:n_H^*","Critical gas density for the formation of a cold interstellar gas phase, \alpha=(n_H/n_H^*)^\beta");  
	control_parameter_add2(control_parameter_double,&agn_parameters.bondi_exponent,"agn:bondi_exponent","agn:beta","exponential dependence of normalization of Bondi-Hoyle accretion rate on ISM gas density"); 
	control_parameter_add2(control_parameter_double,&agn_parameters.minimum_feedback_temperature,"agn:minimum_feedback_temperature","agn:Tmin","Minimum temperature that defines critical energy for a heating event to be triggered, E_crit=m_g*k*Temperature/(mu*Mh*(gamma-1))"); 
	control_parameter_add(control_parameter_int, &agn_feedback_storage, "agn:feedback_storage", "Feedback storage for AGN: 0 = no storage, 1 = storage with E_crit = m_g*k*Tmin/(gamma-1)/mu/mH"); 
	control_parameter_add2(control_parameter_double, &agn_parameters.sink_particle_delta, "agn:sink_particle_delta", "agn:delta","Determines accretion radius for sink particle which is proportional to the resolution, r_{acc}=\\delta \\Delta x"); 
	control_parameter_add2(control_parameter_double, &agn_parameters.agn_merge_radius, "agn:merge_radius", "agn:merge_radius","Determines maximum relative distance for two merging black holes with default, multiples of the resolution size"); 
	control_parameter_add2(control_parameter_double, &agn_parameters.agn_merge_velocity, "agn:merge_velocity", "agn:merge_velocity","Determines maximum relative velocity for two merging black holes"); 
	control_parameter_add(control_parameter_int, &agn_dv_on_off, "agn:dv_on_off", "dv on or off for AGN: 0 = no dv in calculation, 1 = dv in calculation"); 
}

void config_verify_agn() {
	cart_assert( agn_parameters.eddington_factor >= 1.0 );
	cart_assert( agn_parameters.bondi_normalization >= 1.0 );
	cart_assert( agn_parameters.radiative_efficiency > 0.0 && agn_parameters.radiative_efficiency <= 1.0 );
	cart_assert( agn_parameters.feedback_efficiency >= 0.0 && agn_parameters.feedback_efficiency <= 1.0 );
	cart_assert( agn_accretion_recipe == 0 || agn_accretion_recipe == 1 || agn_accretion_recipe == 2); 
	cart_assert( agn_merge_radius_type == 0 || agn_merge_radius_type == 1 );
	cart_assert( agn_parameters.critical_SF_gas_density > 0.0 ); 
	cart_assert( agn_parameters.bondi_exponent > 0.0 ); 
	cart_assert( agn_parameters.minimum_feedback_temperature >= 1e6 ); 
	cart_assert( agn_feedback_storage == 0 || agn_feedback_storage == 1 ); 
	cart_assert( agn_parameters.sink_particle_delta > 0.0); 
	if ( agn_merge_radius_type==0 ) {
	  cart_assert( agn_parameters.agn_merge_radius >= 0.0 && agn_parameters.agn_merge_radius <= pow( 2.0, max_level-1 )); 
	}
	cart_assert( agn_parameters.agn_merge_velocity > 0.0 );
	cart_assert( agn_dv_on_off == 0 || agn_dv_on_off == 1 );
}

double Meddfact;
double Mbondifact;
double Efbfact;
double Eagncritfact; 
double crit_SF_gas_density;

#define SINK_ACCRETION_RADIUS	(3.0) /* in units of accretion kernel scale r_K */
#define SINK_RADIUS_SAMPLES		(3)   /* number of sampling points over accretion radius 
                                       * (this is what we'd change to 0 if we wanted to do single cell accretion)  
                                       * normally set to 3 */  
#define SINK_SAMPLES_1D			(2*SINK_RADIUS_SAMPLES+1)                         /* total number of sampling points in 1d */
#define SINK_SAMPLES_3D			(SINK_SAMPLES_1D*SINK_SAMPLES_1D*SINK_SAMPLES_1D) /* number of sample points per accretion radius */

typedef struct {
	int agn_index;
	int agn_cell;
	float sink_volume;

	int sink_cell_proc[SINK_SAMPLES_3D];
	int sink_cell_index[SINK_SAMPLES_3D];

	float sink_density[SINK_SAMPLES_3D];
	float sink_energy[SINK_SAMPLES_3D];
	float sink_momentum[SINK_SAMPLES_3D][nDim];
} agn_sink_data;

void agn_feedback( int level );
void agn_collect_sink_values( int num_level_agn, agn_sink_data *agn_list , int level );
void agn_update_sink_values( int num_level_agn, agn_sink_data *agn_list , int level );
void agn_compute_agn_physics( int num_level_agn, agn_sink_data *agn_list , int level );
double compute_r_K( int ipart, int icell );
void agn_find_mergers();
void agn_merge_list( int ipart, int *merge_list, int merge_count );
int sort_particles_by_mass( const void *, const void *);

double sink_sample_weights[SINK_SAMPLES_3D];
double sink_sample_offset[SINK_SAMPLES_3D][nDim];
int num_sink_samples;

void init_agn() {
	int i, j, k;
	double x, y, z, r2;
	double weight_sum;

	/* compute sink formalism sampled points' weights and offsets */
	num_sink_samples = 0;
	for ( i = -SINK_RADIUS_SAMPLES; i <= SINK_RADIUS_SAMPLES; i++ ) {
		x = SINK_ACCRETION_RADIUS*(double)i;
		for ( j = -SINK_RADIUS_SAMPLES; j <= SINK_RADIUS_SAMPLES; j++ ) {
			y = SINK_ACCRETION_RADIUS*(double)j;
			for ( k = -SINK_RADIUS_SAMPLES; k <= SINK_RADIUS_SAMPLES; k++ ) {
				z = SINK_ACCRETION_RADIUS*(double)k;

				r2 = x*x+y*y+z*z;
				if ( r2 <= SINK_ACCRETION_RADIUS*SINK_ACCRETION_RADIUS ) {
					sink_sample_offset[num_sink_samples][0] = x;
					sink_sample_offset[num_sink_samples][1] = y;
					sink_sample_offset[num_sink_samples][2] = z;

					sink_sample_weights[num_sink_samples] = exp( -r2 );
					num_sink_samples++;
				}
			}
		}
	}

	weight_sum = 0.0;
	for ( i = 0; i < num_sink_samples; i++ ) {
		weight_sum += sink_sample_weights[i];
	}

	for ( i = 0; i < num_sink_samples; i++ ) {
		sink_sample_weights[i] /= weight_sum;
	}

}

void agn_feedback( int level ) {
	/* primary entry point for agn module */
	int i;
	int ipart;
	int num_level_cells, *level_cells;
	int num_level_agn;
	agn_sink_data *agn_list;
		
	double t_next = tl[level] + dtl[level]; 

	/* setup feedback parameters */
	
	/* compute redshift-dependent factors */
	crit_SF_gas_density = agn_parameters.critical_SF_gas_density/constants->XH/units->number_density;

	Meddfact = agn_parameters.eddington_factor*4.0*M_PI*constants->G*
		constants->mp/agn_parameters.radiative_efficiency/constants->sigmaT/constants->c*units->time; 
	Mbondifact = 4.0*M_PI*0.52*(constants->G/units->length/units->length/units->length*units->mass*units->time*units->time)*
		(constants->G/units->length/units->length/units->length*units->mass*units->time*units->time);
	Efbfact = agn_parameters.radiative_efficiency*agn_parameters.feedback_efficiency*(constants->c/units->velocity)*
		(constants->c/units->velocity)*cell_volume_inverse[level];

	switch ( agn_feedback_storage ) { // agn_feedback_storage = 0 to turn off, 1 to turn on
		case 0:
			Eagncritfact = 0.0; 
			break;
		case 1: 
			Eagncritfact = agn_parameters.minimum_feedback_temperature/units->temperature/constants->wmu/(constants->gamma-1)*cell_volume_inverse[level]; 
			/* Note: mH AND k (kB) is included in units->temperature 
             * Also, to match up in comparison to Efb/cell_mass, need to mult. by cell_volume_inverse[level] (??) */
			break;
		default: 
			cart_error("Invalid agn_feedback_storage value: %d", agn_feedback_storage );
	}

	/* select black hole particles to be moved */
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
	int num_agn_cell_vars = 5;
	double r_K = 0.0; 
	double sample_x[nDim];
	int coords[nDim];
	int num_send_requests = 0;	
	int sync_flag[MAX_PROCS];
	int recv_sample_count[MAX_PROCS];
	int send_sample_count[MAX_PROCS];
	MPI_Request send_requests[4*MAX_PROCS];
	MPI_Request recv_sample_requests[MAX_PROCS];
	MPI_Request recv_cell_values_requests[MAX_PROCS];
	MPI_Request recv_cell_indices_requests[MAX_PROCS];
	int *sample_agn_index[MAX_PROCS];
	int *sample_agn_point[MAX_PROCS];
	double *send_positions[MAX_PROCS];
	double *recv_positions;
	int *send_cell_indices[MAX_PROCS];
	int *recv_cell_indices[MAX_PROCS];
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

			icell = cell_find_position( sample_x );

			if ( icell == -1 || !cell_is_local(icell) ) {
				for ( k = 0; k < nDim; k++ ) {
					coords[k] = (int)sample_x[k];
				}

				proc = processor_owner( sfc_index( coords ) );
				send_sample_count[proc]++;
				agn_list[i].sink_cell_proc[j] = proc;
			} else {
				agn_list[i].sink_cell_proc[j] = local_proc_id;
				agn_list[i].sink_cell_index[j] = icell;

				agn_list[i].sink_density[j] = cell_gas_density(icell);
				agn_list[i].sink_energy[j] = cell_gas_sound_speed(icell); 

				for ( k = 0; k < nDim; k++ ) {
					agn_list[i].sink_momentum[j][k] = cell_momentum(icell, k);
				}
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
			recv_cell_values[i] = cart_alloc(float, num_agn_cell_vars*send_sample_count[i] );

			MPI_Irecv( recv_cell_indices[i], send_sample_count[i], MPI_INT,
					i, 1, mpi.comm.run, &recv_cell_indices_requests[i] );
			MPI_Irecv( recv_cell_values[i], num_agn_cell_vars*send_sample_count[i], 
					MPI_FLOAT, i, 1, mpi.comm.run, 
					&recv_cell_values_requests[i] );
		} else {
			recv_cell_indices_requests[i] = MPI_REQUEST_NULL;
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
				send_cell_values[i] = cart_alloc( float, num_agn_cell_vars*recv_sample_count[i] );

				/* process received sample requests */
				for ( j = 0; j < recv_sample_count[i]; j++ ) {
					icell = cell_find_position( &recv_positions[nDim*j] ); 
					cart_assert( icell != NULL_OCT );

					send_cell_indices[i][j] = icell;
					send_cell_values[i][num_agn_cell_vars*j] = cell_gas_density(icell);
					send_cell_values[i][num_agn_cell_vars*j+1] = cell_gas_sound_speed(icell);
					send_cell_values[i][num_agn_cell_vars*j+2] = cell_momentum(icell,0);
					send_cell_values[i][num_agn_cell_vars*j+3] = cell_momentum(icell,1);
					send_cell_values[i][num_agn_cell_vars*j+4] = cell_momentum(icell,2);
				}

				MPI_Isend( send_cell_indices[i], recv_sample_count[i], MPI_INT,
						i, 1, mpi.comm.run, &send_requests[num_send_requests++] );
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

			for ( j = 0; j < send_sample_count[i]; j++ ) {
				agn_index = sample_agn_index[i][j];
				sample = sample_agn_point[i][j];

				agn_list[agn_index].sink_cell_index[sample] = recv_cell_indices[i][j];

				agn_list[agn_index].sink_density[sample] = recv_cell_values[i][num_agn_cell_vars*j];
				agn_list[agn_index].sink_energy[sample] = recv_cell_values[i][num_agn_cell_vars*j+1];
				agn_list[agn_index].sink_momentum[sample][0] = recv_cell_values[i][num_agn_cell_vars*j+2];
				agn_list[agn_index].sink_momentum[sample][1] = recv_cell_values[i][num_agn_cell_vars*j+3];
				agn_list[agn_index].sink_momentum[sample][2] = recv_cell_values[i][num_agn_cell_vars*j+4];
			}
			cart_free( recv_cell_indices[i] );
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
	MPI_Request send_requests[4*MAX_PROCS];
	MPI_Request recv_sample_requests[MAX_PROCS];
	int *send_cell_indices[MAX_PROCS];
	float *send_cell_deltas[MAX_PROCS];
	int *recv_cell_indices[MAX_PROCS];
	float *recv_cell_deltas[MAX_PROCS];
	int min_level_changed = level;
	int max_level_changed = level;

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
				/* Note that sink_density here has units of mass (replaced list quantities) */
				max_level_changed = max(max_level_changed, cell_level(icell));
				min_level_changed = min(min_level_changed, cell_level(icell));

				density_fraction = new_density / cell_gas_density(icell);

				cell_gas_density(icell) = new_density;
				cell_gas_energy(icell) *= density_fraction;
				cell_gas_internal_energy(icell) *= density_fraction;

				for ( k = 0; k < nDim; k++ ) {
					cell_momentum(icell,k) = density_fraction*cell_momentum(icell,k) + agn_list[i].sink_momentum[j][k];
				}

				/* energy due to feedback */
				cell_gas_energy(icell) += agn_list[i].sink_energy[j];
				cell_gas_internal_energy(icell) += agn_list[i].sink_energy[j];
				cell_gas_pressure(icell) = max( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );
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

					/* energy due to feedback */
					cell_gas_energy(icell) += recv_cell_deltas[i][j*num_agn_cell_vars+1];
					cell_gas_internal_energy(icell) += recv_cell_deltas[i][j*num_agn_cell_vars+1];
					cell_gas_pressure(icell) = max( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );
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

#ifdef DEBUG
	if ( max_level_changed > level || min_level_changed < level ) {
	  cart_debug("level = %d, max_level_changed = %d, min_level_changed = %d, timestep = %e ", level, max_level_changed, min_level_changed, tl[level]);
	}
#endif
}


void agn_compute_agn_physics( int num_level_agn, agn_sink_data *agn_list, int level  ) {
	int i,j,k;
	int ipart;
	int agn_cell;
	double cs_ideal, pressure_ideal;
	double sink_averaged_cs, sink_averaged_density; 
	double sink_averaged_gas_velocity[nDim];
	double change_in_energy[SINK_SAMPLES_3D], change_in_mass[SINK_SAMPLES_3D], change_in_momentum[SINK_SAMPLES_3D][nDim]; 
	double Macc, newMacc, Medd, new_momentum[nDim], new_particle_mass, Mbondi, Mbondibetafact, Mbh, Efb, virtual_sample_volume, sink_mass, dv = 0.0, dv_cell;
	double delta_t;

	double t_next = tl[level] + dtl[level]; 
	double Efb_print;

#ifdef COSMOLOGY 
#pragma omp parallel for default(none), private( cs_ideal, pressure_ideal, i, j, k, ipart, Mbh, virtual_sample_volume, delta_t, sink_mass, sink_averaged_cs, sink_averaged_density, sink_averaged_gas_velocity, new_momentum, Medd, Macc, dv, Mbondibetafact, Mbondi, change_in_mass, newMacc, change_in_momentum, new_particle_mass, Efb, Efb_print, agn_cell, dv_cell ), shared( num_level_agn,  agn_list, t_next, num_sink_samples, sink_sample_weights, Meddfact, Mbondifact, Efbfact, Eagncritfact, agn_accretion_recipe, particle_x, particle_v, particle_t, agn_parameters, units, constants, particle_mass, star_metallicity_II, cell_vars, tl, abox, level, crit_SF_gas_density, agn_dv_on_off, particle_id )
#else
#pragma omp parallel for default(none), private( cs_ideal, pressure_ideal, i, j, k, ipart, Mbh, virtual_sample_volume, delta_t, sink_mass, sink_averaged_cs, sink_averaged_density, sink_averaged_gas_velocity, new_momentum, Medd, Macc, dv, Mbondibetafact, Mbondi, change_in_mass, newMacc, change_in_momentum, new_particle_mass, Efb, Efb_print, agn_cell, dv_cell ), shared( num_level_agn,  agn_list, t_next, num_sink_samples, sink_sample_weights, Meddfact, Mbondifact, Efbfact, Eagncritfact, agn_accretion_recipe, particle_x, particle_v, particle_t, agn_parameters, units, constants, particle_mass, star_metallicity_II, cell_vars, tl, level, crit_SF_gas_density, agn_dv_on_off, particle_id )
#endif
	/* loop over black holes */
	for ( i=0; i<num_level_agn; i++ ) {
		ipart = agn_list[i].agn_index;
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

				Mbondi = agn_parameters.bondi_normalization*Mbondifact * Mbh*Mbh * sink_averaged_density * 
							pow( agn_dv_on_off*dv + sink_averaged_cs*sink_averaged_cs, -1.5 )*delta_t; 
				Macc = min( Medd, Mbondi );
				break;
			case 2: /* Bondi accretion (constant beta) */
				dv = (particle_v[ipart][0]-sink_averaged_gas_velocity[0])*(particle_v[ipart][0]-sink_averaged_gas_velocity[0])+
					(particle_v[ipart][1]-sink_averaged_gas_velocity[1])*(particle_v[ipart][1]-sink_averaged_gas_velocity[1])+
					(particle_v[ipart][2]-sink_averaged_gas_velocity[2])*(particle_v[ipart][2]-sink_averaged_gas_velocity[2]);

				if ( sink_averaged_density > crit_SF_gas_density ) {
					Mbondibetafact = pow( sink_averaged_density / ( crit_SF_gas_density ), agn_parameters.bondi_exponent ) * Mbondifact; 
				} else {
					Mbondibetafact = Mbondifact;
				}
				Mbondi = Mbondibetafact * Mbh*Mbh * sink_averaged_density * pow( agn_dv_on_off*dv + sink_averaged_cs*sink_averaged_cs, -1.5 )*delta_t;
				Macc = min( Medd, Mbondi );
				break;
			default :
				cart_error("Invalid agn_accretion_recipe in agn_compute_agn_physics");
		}

		/* Calculate change_in_energy (due to feedback), change_in_mass, change_in_momentum */

		for ( j=0; j<num_sink_samples; j++ ) { 
			change_in_mass[j] = min( Macc * sink_sample_weights[j] , .667 * agn_list[i].sink_density[j] * virtual_sample_volume ); 
			newMacc += change_in_mass[j];

			for ( k=0; k<nDim; k++ ) {
				change_in_momentum[j][k] = ( change_in_mass[j] / agn_list[i].sink_density[j] ) * agn_list[i].sink_momentum[j][k]; 
				new_momentum[k] += change_in_momentum[j][k];
			}
		}

		/* Calculate new velocity and mass for black hole from accreted gas mass */
		new_particle_mass = particle_mass[ipart] + ( 1-agn_parameters.radiative_efficiency ) * newMacc;

		for ( j=0; j<nDim; j++ ) {
			particle_v[ipart][j] = ( particle_v[ipart][j] * particle_mass[ipart] + new_momentum[j] ) / new_particle_mass;
		}

		particle_mass[ipart] = new_particle_mass;

		/* Compute feedback energy and change_in_energy */

		Efb = Efbfact*newMacc;
		Efb_print = Efb;
#ifndef ENRICH 
#error "Enrichment required for AGN stored feedback energy"
#endif /* ENRICH */
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
    
double compute_r_K( int ipart, int icell ) {  
	double r_K;
	double dv_cell;
	double r_BH; /* Bondi-Hoyle radius */
	double Delta_x=cell_size[max_level];
	double Mbh = particle_mass[ipart];

	cart_assert( ipart >= 0 && ipart < num_star_particles );
	cart_assert( icell >= 0 && icell < num_cells );

	dv_cell =    (particle_v[ipart][0]-cell_momentum(icell,0)/cell_gas_density(icell)) *
		(particle_v[ipart][0]-cell_momentum(icell,0)/cell_gas_density(icell)) +
		(particle_v[ipart][1]-cell_momentum(icell,1)/cell_gas_density(icell)) *
		(particle_v[ipart][1]-cell_momentum(icell,1)/cell_gas_density(icell)) +
		(particle_v[ipart][2]-cell_momentum(icell,2)/cell_gas_density(icell)) *
		(particle_v[ipart][2]-cell_momentum(icell,2)/cell_gas_density(icell));   

	r_BH = constants->G * Mbh / ( cell_gas_sound_speed(icell) * cell_gas_sound_speed(icell) + dv_cell );
	if ( r_BH < Delta_x/4 ) {
		r_K = Delta_x/4;
	}

	if (  (r_BH >= Delta_x/4) && (r_BH <= agn_parameters.sink_particle_delta * Delta_x/2) )  {
		r_K = r_BH;
	}

	if ( r_BH > agn_parameters.sink_particle_delta * Delta_x/2 ) {
		r_K = agn_parameters.sink_particle_delta * Delta_x/2;
	}

	return r_K;
}

void agn_find_mergers() {
	int i, j, k;
	int dx, di, dj, dk;
	int coords[nDim];
	int coords2[nDim];
	int icell, ipart, ipart2;
	int num_agn, num_root;
	int agn_merger_count;
	int *agn_particle_list;
	int *agn_index;
	int *num_root_agn;
	int *order;

	double max_merging_distance;

#define MAX_AGN_MERGERS		10
	int merge_list[MAX_AGN_MERGERS];

	/* create buffer of agn particles */
	build_particle_buffer( num_particle_species - 1, STAR_TYPE_AGN );

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
			cart_assert( icell >= 0 && icell < num_cells_per_level[min_level]+num_buffer_cells[min_level] );
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
	order = cart_alloc( int, num_agn );

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL &&
				particle_is_star(ipart) && star_particle_type[ipart] == STAR_TYPE_AGN ) {
			icell = cell_find_position_level( min_level, particle_x[ipart] );
			cart_assert( icell >= 0 && icell < num_cells_per_level[min_level]+num_buffer_cells[min_level] );
			agn_particle_list[ agn_index[icell] + num_root_agn[icell]++ ] = ipart;
		}
	}

	for ( i = 0; i < num_root-1; i++ ) {
		cart_assert( agn_index[i]+num_root_agn[i] == agn_index[i+1] );
	}

	j = 0;
	for ( i = 0; i < num_root; i++ ) {
		j += num_root_agn[i];
	}
	cart_assert( j == num_agn );	

	/* sort by mass */
	for ( i = 0; i < num_agn; i++ ) {
		order[i] = agn_particle_list[i];
	}

	qsort( order, num_agn, sizeof(int), sort_particles_by_mass );

	for ( i = 0; i < num_agn; i++ ) {
		cart_assert( order[i] >= 0 && order[i] < num_particles );
		cart_assert( particle_level[order[i]] != FREE_PARTICLE_LEVEL );
		cart_assert( particle_is_star(order[i]) && star_particle_type[order[i]] == STAR_TYPE_AGN );
	}

	switch ( agn_merge_radius_type ) {
		case 0 :  /* multiples of the resolution size */
			cart_assert( agn_parameters.agn_merge_radius >= 0.0 && agn_parameters.agn_merge_radius*cell_size[max_level] <= 0.5*cell_size[min_level] ); //CA-10/10/11    
			max_merging_distance = agn_parameters.agn_merge_radius * cell_size[max_level];
			break;
		case 1 :  /* comoving Mpc/h */
#ifdef COSMOLOGY                                                                                                                                
           cart_assert( agn_parameters.agn_merge_radius >= 0.0 &&
                    agn_parameters.agn_merge_radius/units->length_in_chimps <= 0.5*cell_size[min_level] );
           max_merging_distance = agn_parameters.agn_merge_radius / units -> length_in_chimps;
#else
           max_merging_distance = agn_parameters.agn_merge_radius * constants->Mpc / units->length;
#endif
		   break;
	}

	dx = (int)(max_merging_distance / cell_size[min_level]) + 1;

	cart_debug("max_merging_distance = %e, dx = %d", max_merging_distance, dx );

	for ( i = num_agn-1; i >= 0; i-- ) {
		ipart = order[i];

		/* skip over AGN that have been removed due to merging with more massive AGN */
		if ( particle_level[ipart] != FREE_PARTICLE_LEVEL ) {
			cart_assert( ipart >= 0 && ipart < num_particles );
			cart_assert( particle_is_star(ipart) && star_particle_type[ipart] == STAR_TYPE_AGN );
			cart_debug("finding mergers for agn %d, id = %d, Mbh = %e Msun", i, particle_id[ipart], particle_mass[ipart] * units->mass / constants->Msun );

			agn_merger_count = 0;

			for ( j = 0; j < nDim; j++ ) {
				coords[j] = (_nt)(particle_x[ipart][j]);
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

								cart_assert( particle_is_star(ipart2) && star_particle_type[ipart2] == STAR_TYPE_AGN );

								/* compare the two black holes */
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
	double particle_mass_initial = particle_mass[ipart];
	double particle_x_initial[nDim];
	double particle_v_initial[nDim];

	for ( j = 0; j < nDim; j++ ) {
		particle_x_initial[j] = particle_x[ipart][j];
		particle_v_initial[j] = particle_v[ipart][j];
	}

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
	int i;
	double max_merging_distance;  
	double net_distance;
	double net_velocity;

	switch ( agn_merge_radius_type ) {
		case 0 :  /* multiples of the resolution size */
			cart_assert( agn_parameters.agn_merge_radius >= 0.0 && 
					agn_parameters.agn_merge_radius*cell_size[max_level] <= 0.5*cell_size[min_level] ); 
			max_merging_distance = agn_parameters.agn_merge_radius * cell_size[max_level];
			break;
		case 1 :  /* comoving Mpc/h */
#ifdef COSMOLOGY
			cart_assert( agn_parameters.agn_merge_radius >= 0.0 && 
					agn_parameters.agn_merge_radius/units->length_in_chimps <= 0.5*cell_size[min_level] );
			max_merging_distance = agn_parameters.agn_merge_radius / units->length_in_chimps;
#else
			max_merging_distance = agn_parameters.agn_merge_radius * constants->Mpc / units->length;
#endif
			break;
	}

	/* compare the relative distance of the two AGN */
	net_distance = compute_distance_periodic( particle_x[ipart], particle_x[ipart2] );

	net_velocity = sqrt( (particle_v[ipart][0]-particle_v[ipart2][0]) *
			(particle_v[ipart][0]-particle_v[ipart2][0]) + 
			(particle_v[ipart][1]-particle_v[ipart2][1]) *
			(particle_v[ipart][1]-particle_v[ipart2][1]) + 
			(particle_v[ipart][2]-particle_v[ipart2][2]) *
			(particle_v[ipart][2]-particle_v[ipart2][2]) );

	if ( net_distance < max_merging_distance && 
			net_velocity*units->velocity/constants->kms < agn_parameters.agn_merge_velocity ) {
		return 1;
	} else {
		return 0;
	}
}

int sort_particles_by_mass( const void *a, const void *b ) {
	int index_a = *(int *)a;
	int index_b = *(int *)b;

	cart_assert( index_a >= 0 && index_a < num_particles );
	cart_assert( index_b >= 0 && index_b < num_particles );

	if ( particle_mass[index_a] < particle_mass[index_b] ) {
		return -1;
	} else if ( particle_mass[index_a] > particle_mass[index_b] ) {
		return 1;
	} else {
		/* decreasing order of id */
		return ( particle_id[index_b] - particle_id[index_a] );
	}
} 

#endif /* STARFORM && AGN */
