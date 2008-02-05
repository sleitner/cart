#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "timestep.h"
#include "tree.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "cell_buffer.h"
#include "parallel.h"
#include "refinement.h"
#include "iterators.h"
#include "viewdump.h"
#include "density.h"
#include "timing.h"
#include "gravity.h"
#include "units.h"
#include "particle.h"
#include "starformation.h"
#include "io.h"
#include "auxiliary.h"
#include "cooling.h"

double tl[max_level-min_level+1];
double tl_old[max_level-min_level+1];
double dtl[max_level-min_level+1];
double dtl_old[max_level-min_level+1];
double aexp[max_level-min_level+1];
double aexp_old[max_level-min_level+1];

int num_steps_on_level[max_level-min_level+1];

double a_init;
double a_end = 1.0;
double t_init;
double t_end = 0.0;

int output_frequency = 0;
int restart_frequency = 1;
int particle_output_frequency = 0;
int tracer_output_frequency = 0;
int grid_output_frequency = 0;
int max_steps = 0;
int max_cfl_sync_level = 2;

double cfl = 0.6;
double particle_cfl = 0.0;
double max_time_inc = 1.2;
double min_time_dec = 1.25;
double max_da = 3e-3;
double max_dt = 0.125;
double max_frac_da = 0.1;

int step;
int step_of_last_increase = 0;
int steps_before_increasing = 4;

int global_timestep( double dt ) {
	int i;
	int ret;
	int global_ret;
	int sf;
	double dtratio;

	/* set old vars */
	for ( i = min_level; i <= max_level; i++ ) {
		tl_old[i] = tl[i];
		dtl_old[i] = dtl[i];
		aexp_old[i] = aexp[i];
		num_steps_on_level[i] = 0;
	}

	dtl[min_level] = dt;

	/* ensure consistency of time variables */
	for ( i = min_level+1; i <= max_level; i++ ) {
		tl[i] = tl[min_level];
		aexp[i] = aexp[min_level];
		dtl[i] = 0.5*dtl[i-1];
	}

#ifdef COOLING
	/* prepare for cooling timestep */
	set_cooling_redshift( aexp[min_level] );
#endif /* COOLING */

#ifdef STARFORM
	num_new_stars = 0;
	last_star_id = particle_species_indices[num_particle_species]-1;

	/* compute frequency of star formation calls */
	for ( i = min_level; i <= max_level; i++ ) {
		dtratio = max( dtmin_SF / ( t0 * aexp[i] * aexp[i] * dtl[i] ), 1e-30 );
		sf = max( 0, nearest_int( log(dtratio)/log(2.0) ) );
		sf = min( sf, i );
		star_formation_frequency[i] = 1 << sf;
	}
#endif /* STARFORM */

	ret = timestep( min_level );

	/* check if any other processors had problems (violation of CFL condition) */
	MPI_Allreduce( &ret, &global_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );

	/* do a last refinement step (without allowing derefinement */
	if ( global_ret != -1 ) {
#ifdef GRAVITY
		for ( i = min_level; i <= max_level-1; i++ ) {
			assign_density(i);
		}
#endif

		for ( i = min_level; i <= max_level-1; i++ ) {
			modify( i, 0 );
		}
#ifdef STARFORM
		/* now remap ids of stars created in this timestep */
		remap_star_ids();
#endif /* STARFORM */
	}

	return global_ret;
}

int timestep( int level ) 
/* returns -1 if timestep would invalidate cfl condition */
{
	int courant_cell;
	double velocity;
	int ret;
	int step;
	int step_ret;
	int true_ret;
	int nlevel;
	double dt_needed;

	cart_assert( level >= min_level && level <= max_level );

	/* assume step was sucessful */
	ret = 0;

	start_timing_level( level );
	start_time( LEVEL_TIMER );

        if ( local_proc_id == MASTER_NODE  ) {
                cart_debug("timestep(%u, %e)", level, dtl[level] ); 
        }


#ifdef HYDRO
	hydro_copy_vars( level, COPY_ZERO_REF, COPY_SPLIT_NEIGHBORS );	
#endif /* HYDRO */

#ifdef GRAVITY
#ifdef PARTICLES
	compute_accelerations_particles(level);
#endif /* PARTICLES */
#endif /* GRAVITY */
	
	if ( level < max_level && level < max_level_now() ) {
		step_ret = timestep( level + 1 );
		ret = min( ret, step_ret );
		if ( ret == -1 && level < max_cfl_sync_level ) { 
			end_time( LEVEL_TIMER );
			end_timing_level( level );
			return ret; 
		}
		step_ret = timestep( level + 1);
		ret = min( ret, step_ret );
		if ( ret == -1 && level < max_cfl_sync_level ) { 
			end_time( LEVEL_TIMER );
			end_timing_level( level );
			return ret; 
		}
	} else {
		/* advance timestep on lower levels */
		for ( nlevel = level + 1; nlevel <= max_level; nlevel++ ) {
			for ( step = 0; step < 1<<(nlevel - level); step++ ) {

				if ( nlevel <= max_cfl_sync_level ) {
					MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );

					if ( true_ret < 0 ) {
						end_time( LEVEL_TIMER );
						end_timing_level( level );
						return true_ret;
					}
				}

				tl_old[nlevel] = tl[nlevel];
				tl[nlevel] += dtl[nlevel];

#ifdef COSMOLOGY
			        aexp_old[nlevel] = aexp[nlevel];
			        aexp[nlevel] = b2a( tl[nlevel] );
#endif

				num_steps_on_level[nlevel]++;
			}

			cart_assert( fabs( tl[nlevel] - (tl[level]+dtl[level]) ) < 1e-6 );
		}
	}

#ifdef GRAVITY
	if ( level > min_level && level < max_level ) {
		restrict_to_level( level );
	}
#endif

#ifdef HYDRO

#ifdef GRAVITY
	interpolate_potential( level );
	compute_accelerations_hydro( level );
#endif /* GRAVITY */
	
	/* test if timestep is still valid */
	hydro_timestep( level, &courant_cell, &velocity );
	dt_needed = cfl * cell_size[level] / velocity;

	/* check for cfl condition violation... */
	if ( (dtl[level]-dt_needed) > 1e-12 ) {
		cart_debug("uh oh, violated cfl condition current dt: %0.25e needed %0.25e courant_cell = %u, velocity = %e", 
			dtl[level], dt_needed, courant_cell, velocity );

		cart_debug("density = %e, pressure = %e, as = %e", cell_gas_density(courant_cell),
			cell_gas_pressure(courant_cell), 
			sqrt( gamma * cell_gas_pressure(courant_cell) / cell_gas_gamma(courant_cell ) ) );
		cart_debug("momentum = %e %e %e", cell_momentum(courant_cell,0), cell_momentum(courant_cell,1),
			cell_momentum(courant_cell,2) );

		ret = -1;
	}

	if ( level <= max_cfl_sync_level ) {
		MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );

		if ( true_ret < 0 ) {
			end_time( LEVEL_TIMER );
			end_timing_level( level ); 
			return true_ret;
		}
	}

	/* do hydro step */
	hydro_step( level );

#endif /* HYDRO */

#ifdef HYDRO_TRACERS
	move_hydro_tracers( level );
	update_tracer_list( level );
#endif /* HYDRO_TRACERS */

#ifdef PARTICLES

#ifdef HYDRO
#ifdef STARFORM
        if ( num_steps_on_level[level] % star_formation_frequency[level] == 0 ) {
                star_formation( level, star_formation_frequency[level] );
        }
#endif /* STARFORM */
#endif /* HYDRO */

#ifdef GRAVITY
#ifdef HYDRO
	/* if hydro is enabled, we destroyed the value of the acceleration variable
	 * and need to recompute it here */
	compute_accelerations_particles(level);
#endif /* HYDRO */
#endif /* GRAVITY */

	move_particles( level );
	update_particle_list( level );

#ifdef STARFORM
	/* update cell values changed by starformation and feedback */
	update_buffer_level( level, all_hydro_vars, num_hydro_vars );
#endif /* STARFORM */

#endif /* PARTICLES */

	/* advance time on level */
	tl_old[level] = tl[level];
	tl[level] += dtl[level];

#ifdef COSMOLOGY
	aexp_old[level] = aexp[level];
	aexp[level] = b2a( tl[level] );
#endif

#ifdef GRAVITY
	/* recompute potential */
	assign_density( level );

#ifdef HYDRO
	copy_potential( level );
#endif

	solve_poisson( level, num_steps_on_level[level] );
#endif /* GRAVITY */

	if ( level < max_level ) {
		modify( level, 1 );
	}

	num_steps_on_level[level]++;

	end_time( LEVEL_TIMER );
	end_timing_level( level );

	return ret;
}

void choose_timestep( double *dt ) {
	int i, j;
	double velocity;
	double level_velocity;
	int courant_cell, level_courant_cell;
        double adum1, adum2, dda;
	double dt_new, dt_min;

#ifdef CONSTANT_TIMESTEP
	return;
#else 

	dt_new = 0.0;

#ifdef HYDRO 
	velocity = 0.0;
	for ( i = min_level; i <= max_level; i++ ) {
		hydro_timestep( i, &level_courant_cell, &level_velocity );
		if ( level_velocity > velocity ) {
			velocity = level_velocity;
			courant_cell = level_courant_cell;
		}
	}

	cart_debug("cfl cell: velocity = %e, cell = %u, level = %u, pressure = %e, density = %e, momentum = %e %e %e",
		velocity, courant_cell, cell_level(courant_cell), cell_gas_pressure(courant_cell),
		cell_gas_density(courant_cell), cell_momentum(courant_cell,0),
		cell_momentum(courant_cell,1), cell_momentum(courant_cell,2) );

	cart_assert( velocity > 0.0 );

	if ( dt_new > 0.0 ) {
		dt_new = min( dt_new, cfl *cell_size[min_level] / velocity );
	} else {
		dt_new = cfl *cell_size[min_level] / velocity;
	}
#endif /* HYDRO */

#ifdef PARTICLES
	if ( particle_cfl != 0.0 ) {
		for ( i = 0; i < num_particles; i++ ) {
			if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
				velocity = 0.0;
				for ( j = 0; j < nDim; j++ ) {
					velocity = max( fabs(particle_v[i][j]), velocity );
				}

				if ( dt_new > 0.0 ) {
					dt_new = min( dt_new, particle_cfl * cell_size[min_level] / velocity );
				} else {
					dt_new = particle_cfl * cell_size[min_level] / velocity;
				}
			}	
		}
	}
#endif /* PARTICLES */

	start_time( COMMUNICATION_TIMER );
	MPI_Reduce( &dt_new, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD );
	end_time( COMMUNICATION_TIMER );

	if ( local_proc_id == MASTER_NODE ) {
		/* if a cfl condition forces a decrease in the timestep,
		 * force it to be at least min_time_dec smaller so we avoid
		 * having to restart with a smaller timestep */
		if ( *dt > 0.0 && (dt_min-(*dt))/(*dt) < 1e-6 ) {
			dt_new = min( (*dt) / min_time_dec, dt_min );
			step_of_last_increase = step;
		} else {
			dt_new = dt_min;
		}

#ifdef COSMOLOGY
	        adum1 = b2a ( tl[min_level] );

		if ( dt_new > 0.0 ) {
			adum2 = b2a ( min( tl[min_level] + dt_new, 0.0 )  );
			dda = min( (adum2 - adum1)/adum1 , max_frac_da );
			
			if ( max_da > 0.0 ) {
				dda = min( dda*adum1 , max_da );
			}

			dt_new = min( dt_new, a2b(adum1+dda)-tl[min_level] );
		} else {
			dt_new = a2b(adum1+max_da)-tl[min_level];
		}
#endif /* COSMOLOGY */

		if ( max_dt > 0.0 ) {
			if ( dt_new > 0.0 ) {
				dt_new = min( dt_new, max_dt );
			} else {
				dt_new = max_dt;
			}
		}

		/* enforce maximum change in timestep (dt=0 for first step) */
		if ( *dt > 0.0 && (dt_new-(*dt))/(*dt) > -1e-6 ) {
			if ( step - step_of_last_increase >= steps_before_increasing ) {
				dt_new = min( (*dt)*max_time_inc, dt_new );
				step_of_last_increase = step;
			} else {
				dt_new = min( (*dt), dt_new );
				if ( step_of_last_increase == 0 ) {
					step_of_last_increase = step;
				}
			}
		}

		cart_debug("chose %e as our next timestep", dt_new );
		*dt = dt_new;
	}

	start_time( COMMUNICATION_TIMER );
	MPI_Bcast( dt, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	end_time( COMMUNICATION_TIMER );

#endif /* CONSTANT_TIMESTEP */
}

#ifdef HYDRO
void hydro_timestep( int level, int *courant_cell, double *velocity ) {
	int i, j;
	int icell;
	int num_level_cells;
	int *level_cells;
	int ivas;
	double vas;
	double rho_r, as, vel;

	vas = 0.0;
	ivas = 0;

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf(icell) ) {
			rho_r = 1.0/cell_gas_density(icell);
			as = sqrt( gamma * rho_r * cell_gas_pressure(icell) );
			vel = min_courant_velocity;
			for ( j = 0; j < nDim; j++ ) {
				if ( fabs(cell_momentum(icell,j)) > vel ) {
					vel = fabs(cell_momentum(icell,j));
				}
			}
			vel = vel*rho_r + as;

			if ( vel >= vas ) {
				ivas = icell;
				vas = vel;
			}
		}
	}
	cart_free( level_cells );
	
	*courant_cell = ivas;
	*velocity = vas;
}
#endif /* HYDRO */
