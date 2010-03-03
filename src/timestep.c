#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cooling.h"
#include "cosmology.h"
#include "density.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "io.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_solver.h"
#include "starformation.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


double t_init = 0.0;
double t_end = 0.0;
DEFINE_LEVEL_ARRAY(double,dtl);
DEFINE_LEVEL_ARRAY(double,dtl_old);
DEFINE_LEVEL_ARRAY(double,tl);
DEFINE_LEVEL_ARRAY(double,tl_old);

#ifdef COSMOLOGY
double auni_init = 1.0e-3;
double auni_end = 1.0;
DEFINE_LEVEL_ARRAY(double,abox);
DEFINE_LEVEL_ARRAY(double,abox_old);
DEFINE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

DEFINE_LEVEL_ARRAY(int,num_steps_on_level);

int max_steps = 0;
double timelimit = 0.0;

int max_cfl_sync_level = min_level;

#ifdef HYDRO 
double cfl_run = 0.6;
double cfl_max = 0.6;     /* max allowed, re-do the timestep if that number is exceeded. */
#endif /* HYDRO */

#ifdef PARTICLES
double particle_cfl = 0.0;
#endif /* PARTICLES */

double max_time_inc = 1.1;
double min_time_dec = 1.25;
double max_dt = 0.125;
#ifdef COSMOLOGY
double max_a_inc = 1.1;
double max_da = 3e-3;
#endif /* COSMOLOGY */


int step;
int step_of_last_increase = 0;
int steps_before_increasing = 4;
int current_step_level = -1;

double min_courant_velocity = 1.0e-6;


/*
// NG: it is not clear how make all these 4 parameters consistent.
// For now keep them completely independent
//
#ifdef COSMOLOGY
void control_parameter_set_aini(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  t_init = tcode_from_auni(auni_init);
}


void control_parameter_set_aend(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  t_end = tcode_from_auni(auni_end);
}
#endif


void control_parameter_set_tini(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
#ifdef COSMOLOGY
  auni_init = auni_from_tcode(t_init);
#endif
}


void control_parameter_set_tend(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
#ifdef COSMOLOGY
  auni_end = auni_from_tcode(t_end);
#endif
}
*/


#ifdef HYDRO 
void control_parameter_set_cflrun(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  if(cfl_max < cfl_run) cfl_max = cfl_run;
}


void control_parameter_set_cflmax(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  if(cfl_max < cfl_run) cfl_max = cfl_run;
}
#endif /* HYDRO  */


#ifdef COSMOLOGY
void control_parameter_set_a_inc(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  if(ind == 2)
    {
      max_a_inc += 1.0; /* backward compatibility */
    }
}
#endif /* COSMOLOGY */


void config_init_timestep()
{
#ifdef COSMOLOGY
  ControlParameterOps control_parameter_aini = { control_parameter_set_double, control_parameter_list_double };
  ControlParameterOps control_parameter_aend = { control_parameter_set_double, control_parameter_list_double };
  ControlParameterOps control_parameter_a_inc = { control_parameter_set_a_inc, control_parameter_list_double };
#endif /* COSMOLOGY */
  ControlParameterOps control_parameter_tini = { control_parameter_set_double, control_parameter_list_double };
  ControlParameterOps control_parameter_tend = { control_parameter_set_double, control_parameter_list_double };
#ifdef HYDRO 
  ControlParameterOps control_parameter_cflrun = { control_parameter_set_cflrun, control_parameter_list_double };
  ControlParameterOps control_parameter_cflmax = { control_parameter_set_cflmax, control_parameter_list_double };
#endif /* HYDRO  */

#ifdef COSMOLOGY
  control_parameter_add3(control_parameter_aini,&auni_init,"auni-start","auni_init","a_init","starting value for the cosmic scale factor.");

  control_parameter_add3(control_parameter_aend,&auni_end,"auni-stop","auni_end","a_end","last value for the cosmic scale factor. The simulation stops if this value is reached.");

  control_parameter_add3(control_parameter_a_inc,&max_a_inc,"max-a-increment","max_a_inc","max_frac_da","the largest factor by which the cosmic scale factor is allowed to increase in one time-step.");

  control_parameter_add2(control_parameter_double,&max_da,"max-da","max_da","maximum allowed step in the cosmic scale factor.");
#endif /* COSMOLOGY */

  control_parameter_add2(control_parameter_tini,&t_init,"time-start","t_init","starting value for the code time variable (in code units).");

  control_parameter_add2(control_parameter_tend,&t_end,"time-stop","t_end","last value for the code time variable (in code units). The simulation stops if this value is reached.");

  control_parameter_add2(control_parameter_int,&max_steps,"num-steps","max_steps","number of time-steps to make. Zero value disables this limit.");

  control_parameter_add2(control_parameter_double,&timelimit,"walltime-limit","timelimit","the limit for the wall-clock time for the current job. Zero value disables this limit.");

  control_parameter_add2(control_parameter_int,&max_cfl_sync_level,"max-cfl-sync-level","max_cfl_sync_level","maximum level at which the CFL conditions are synchronized across separate nodes.");

#ifdef HYDRO 
  control_parameter_add3(control_parameter_cflrun,&cfl_run,"cfl-run","cfl_run","cfl","the CFL number for setting the time-step.");

  control_parameter_add3(control_parameter_cflmax,&cfl_max,"cfl-max","cfl_max","cfl","the maximum acceptable CFL number. If this number is exceeded, the time-step needs to be redone. It is a good sense to set this number just a little bit higher than <CFL-run>, to avoid extra restarts due to numerical noise.");
#endif /* HYDRO */

#ifdef PARTICLES
  control_parameter_add2(control_parameter_double,&particle_cfl,"particle-cfl","particle_cfl","the CFL number for particle dynamics. In HYDRO mode this number is usually not needed, as the grid CFL conditions superceeds that of particles. Setting it to zero disables this limit.");
#endif /* PARTICLES */

  control_parameter_add2(control_parameter_double,&max_time_inc,"max-timestep-increment","max_time_inc","the largest factor by which the time-step is allowed to increase.");

  control_parameter_add2(control_parameter_double,&min_time_dec,"min-timestep-decrement","min_time_dec","the smallest factor by which the time-step is allowed to decrease.");

  control_parameter_add2(control_parameter_double,&max_dt,"max-dt","max_dt","maximum allowed time-step.");
}


void config_verify_timestep()
{
#ifdef COSMOLOGY
  cart_assert(auni_init>0.0 && !(auni_init>auni_end));

  cart_assert(max_a_inc > 1.0);

  cart_assert(max_da >= 0.0);
#endif /* COSMOLOGY */

  cart_assert(!(t_init > t_end));

  cart_assert(max_steps >= 0);

  cart_assert(max_cfl_sync_level >= min_level);

#ifdef HYDRO 
  cart_assert(cfl_run>0.0 && !(cfl_run>cfl_max));
#endif /* HYDRO */

#ifdef PARTICLES
  cart_assert(!(particle_cfl < 0.0));
#endif /* PARTICLES */

  cart_assert(max_time_inc > 1.0);

  cart_assert(min_time_dec > 1.0);

  cart_assert(max_dt >= 0.0);
}


int global_timestep( double dt ) {
	int i;
	int ret;
	int global_ret;
	int sf;
	double dtratio;

	start_time( LEVEL_TIMER );
	start_time( WORK_TIMER );

	/* set old vars */
	for ( i = min_level; i <= max_level; i++ ) {
		tl_old[i] = tl[i];
		dtl_old[i] = dtl[i];
#ifdef COSMOLOGY
		abox_old[i] = abox[i];
#endif /* COSMOLOGY */
		num_steps_on_level[i] = 0;
	}

	dtl[min_level] = dt;

	/* ensure consistency of time variables */
	for ( i = min_level+1; i <= max_level; i++ ) {
		tl[i] = tl[min_level];
		dtl[i] = 0.5*dtl[i-1];
#ifdef COSMOLOGY
		abox[i] = abox[min_level];
		auni[i] = auni[min_level];
#endif /* COSMOLOGY */
	}

	units_update(min_level);

	end_time( WORK_TIMER );

#ifdef RADIATIVE_TRANSFER
	rtStepBegin();
	/* By default update tables once per step */
	rtUpdateTables();
#else
#ifdef COOLING
#ifdef COSMOLOGY
	/* prepare for cooling timestep */
	set_cooling_redshift( auni[min_level] );
#else
#error "Setting COOLING and !RADIATIVE_TRANSFER requires COSMOLOGY."
#endif /* COSMOLOGY */
#endif /* COOLING */
#endif /* RADIATIVE_TRANSFER */

#ifdef STARFORM

	start_time( WORK_TIMER );

	num_new_stars = 0;
	last_star_id = particle_species_indices[num_particle_species]-1;

	/* compute frequency of star formation calls */
	for ( i = min_level; i <= max_level; i++ ) {
		dtratio = max( sf_sampling_timescale * constants->yr / units->time , 1e-30 );
		sf = max( 0, nearest_int( log(dtratio)/log(2.0) ) );
		sf = min( sf, i );
		star_formation_frequency[i] = 1 << sf;
	}

	end_time( WORK_TIMER );

#endif /* STARFORM */

	ret = timestep( min_level, MPI_COMM_WORLD );
	current_step_level = -1;

	/* check if any other processors had problems (violation of CFL condition) */
	start_time( COMMUNICATION_TIMER );
	MPI_Allreduce( &ret, &global_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
	end_time( COMMUNICATION_TIMER );

	/* do a last refinement step (without allowing derefinement */
	if ( global_ret != -1 ) {
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		for ( i = min_level; i <= max_level-1; i++ ) {
			assign_density(i);
		}
#endif

#ifdef REFINEMENT
		for ( i = min_level; i <= max_level-1; i++ ) {
			modify( i, 0 );
		}
#endif

#ifdef STARFORM
		/* now remap ids of stars created in this timestep */
		remap_star_ids();
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER	
		rtStepEnd();
#endif /* RADIATIVE_TRANSFER */
	}

	end_time(LEVEL_TIMER);
	return global_ret;
}

int timestep( int level, MPI_Comm local_comm ) 
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
	MPI_Comm child_comm;
	int refined;

	cart_assert( level >= min_level && level <= max_level );

	current_step_level = level;

	/* assume step was sucessful */
	ret = 0;

	start_timing_level( level );
	start_time( LEVEL_TIMER );
	
	start_time( COMMUNICATION_TIMER );
        if(mpi_custom_flags & MPI_CUSTOM_SYNC)
          {
            /*
            //  Sync all procs at this level before creating a new communicator
            */
            MPI_Barrier(local_comm);
          }

	/* 
	//  Create a child communicator
	*/
	refined = (level<max_level && level<max_level_now());
	MPI_Comm_split(local_comm,refined,local_proc_id,&child_comm);

	end_time( COMMUNICATION_TIMER );

#ifdef HYDRO
	hydro_copy_vars( level, COPY_ZERO_REF, COPY_SPLIT_NEIGHBORS );	
#endif /* HYDRO */

#if defined(GRAVITY) && defined(PARTICLES)
	compute_accelerations_particles(level);
#endif /* GRAVITY && PARTICLES */
	
	start_time( LOWER_LEVEL_TIMER );  /* this is for internal accounting only */

	if ( level < max_level && level < max_level_now() ) {
		step_ret = timestep( level + 1, child_comm );
		current_step_level = level;
		ret = min( ret, step_ret );
		if ( ret == -1 && level < max_cfl_sync_level ) { 
			end_time( LOWER_LEVEL_TIMER );
			end_time( LEVEL_TIMER );
			end_timing_level( level );
			MPI_Comm_free( &child_comm );
			return ret; 
		}
		step_ret = timestep( level + 1, child_comm );
		current_step_level = level;
		ret = min( ret, step_ret );
		if ( ret == -1 && level < max_cfl_sync_level ) { 
			end_time( LOWER_LEVEL_TIMER );
			end_time( LEVEL_TIMER );
			end_timing_level( level );
			MPI_Comm_free( &child_comm );
			return ret; 
		}
	} else {
		/* advance timestep on lower levels */
		for ( nlevel = level + 1; nlevel <= max_level; nlevel++ ) {
			for ( step = 0; step < 1<<(nlevel - level); step++ ) {

				if ( nlevel <= max_cfl_sync_level ) {
					start_time( COMMUNICATION_TIMER );
					MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
					end_time( COMMUNICATION_TIMER );

					if ( true_ret < 0 ) {
						end_time( LEVEL_TIMER );
						end_timing_level( level );
						MPI_Comm_free( &child_comm );
						return true_ret;
					}
				}

				tl_old[nlevel] = tl[nlevel];
				tl[nlevel] += dtl[nlevel];

#ifdef COSMOLOGY
				abox_old[nlevel] = abox[nlevel];
				abox[nlevel] = abox_from_tcode( tl[nlevel] );
				auni[nlevel] = auni_from_tcode( tl[nlevel] );
#endif

				num_steps_on_level[nlevel]++;
			}

			cart_assert( fabs( tl[nlevel] - (tl[level]+dtl[level]) ) < 1e-6 );
		}
	}

	end_time( LOWER_LEVEL_TIMER );

	start_time( WORK_TIMER );
	units_update(level);
	end_time( WORK_TIMER );

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
	start_time( WORK_TIMER );
	hydro_timestep( level, &courant_cell, &velocity );
	dt_needed = cfl_max * cell_size[level] / velocity;

	/* check for cfl condition violation... */
	if ( dtl[level] > dt_needed && ret != -1 ) {
		cart_debug("uh oh, violated cfl condition current dt: %0.25e needed %0.25e courant_cell = %u, velocity = %e", 
			dtl[level], dt_needed, courant_cell, velocity );

		cart_debug("density = %e, pressure = %e, as = %e", cell_gas_density(courant_cell),
			cell_gas_pressure(courant_cell), 
			sqrt( cell_gas_gamma(courant_cell ) * cell_gas_pressure(courant_cell) / cell_gas_density(courant_cell ) ) );
		cart_debug("momentum = %e %e %e", cell_momentum(courant_cell,0), cell_momentum(courant_cell,1),
			cell_momentum(courant_cell,2) );

		ret = -1;
	}
	end_time( WORK_TIMER );

	if ( level <= max_cfl_sync_level ) {
		start_time( COMMUNICATION_TIMER );
		MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
		end_time( COMMUNICATION_TIMER );

		if ( true_ret < 0 ) {
			end_time( LEVEL_TIMER );
			end_timing_level( level ); 
			MPI_Comm_free( &child_comm );
			return true_ret;
		}
	}

#ifdef RADIATIVE_TRANSFER
	/* Do RT step */
	rtLevelUpdate(level,local_comm);
#endif /* RADIATIVE_TRANSFER */

	/* do hydro step */
	hydro_step( level, local_comm );

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
#if defined(HYDRO) || defined(REFINEMENT)
	/* if hydro or refinement are enabled, we destroyed the value of the 
	 * acceleration variable and need to recompute it here */
	compute_accelerations_particles(level);
#endif /* HYDRO */
#endif /* GRAVITY */

	move_particles( level );
	update_particle_list( level );

#ifdef STARFORM
	/* update cell values changed by starformation and feedback */
	start_time( STELLAR_FEEDBACK_UPDATE_TIMER );
	update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	end_time( STELLAR_FEEDBACK_UPDATE_TIMER );
#endif /* STARFORM */

#endif /* PARTICLES */

	/* advance time on level */
	tl_old[level] = tl[level];
	tl[level] += dtl[level];

#ifdef COSMOLOGY
	abox_old[level] = abox[level];
	abox[level] = abox_from_tcode( tl[level] );
	auni[level] = auni_from_tcode( tl[level] );
#endif

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
	assign_density( level );
#endif

#ifdef GRAVITY
	/* recompute potential */

#ifdef HYDRO
	copy_potential( level );
#endif

	solve_poisson( level, num_steps_on_level[level] );
#endif /* GRAVITY */

#ifdef REFINEMENT
	if ( level < max_level ) {
		modify( level, 1 );
	}
#endif

	cart_debug("timestep(%u, %e, %d)", level, dtl[level], num_steps_on_level[level] );

	num_steps_on_level[level]++;

	end_time( LEVEL_TIMER );
	end_timing_level( level );

	MPI_Comm_free( &child_comm );

	return ret;
}


void choose_timestep( double *dt ) {
	int i, j;
	double velocity;
	double level_velocity;
	int courant_cell, level_courant_cell;
	double adum1, dda;
	double dt_new, dt_min;

	start_time( CHOOSE_TIMESTEP_TIMER );

#ifdef CONSTANT_TIMESTEP

	if ( *dt == 0.0 ) {
		*dt = max_dt;
	}

#ifdef RADIATIVE_TRANSFER
	rtModifyTimeStep(dt);
#endif /* RADIATIVE_TRANSFER */

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
		dt_new = min( dt_new, cfl_run *cell_size[min_level] / velocity );
	} else {
		dt_new = cfl_run *cell_size[min_level] / velocity;
	}
#endif /* HYDRO */

#ifdef PARTICLES
	if ( particle_cfl > 0.0 ) {
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

	start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
	MPI_Reduce( &dt_new, &dt_min, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD );
	end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );

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
		adum1 = abox_from_tcode( tl[min_level] );

		if ( dt_new > 0.0 ) {
			dda = abox_from_tcode( min( tl[min_level] + dt_new, tcode_from_abox(adum1*max_a_inc) ) ) - adum1;  /* dt_new could be very large, it can break cosmology module */
			if ( max_da > 0.0 ) {
				dda = min( dda , max_da );
			}

			dt_new = min( dt_new, tcode_from_abox(adum1+dda)-tl[min_level] );
		} else {
			dt_new = tcode_from_abox(adum1+max_da) - tl[min_level];
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

	start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
	MPI_Bcast( dt, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );	

#endif /* CONSTANT_TIMESTEP */

	end_time( CHOOSE_TIMESTEP_TIMER );
}

#ifdef HYDRO
void hydro_timestep( int level, int *courant_cell, double *velocity ) {
	int i, j;
	int icell;
	int num_level_cells;
	int *level_cells;
	int ivas;
	double vas;
	double as, vel;

	vas = 0.0;
	ivas = 0;

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf(icell) ) {
			as = sqrt(cell_gas_gamma(icell)*cell_gas_pressure(icell)/cell_gas_density(icell));
			vel = cell_gas_density(icell)*min_courant_velocity;
			for ( j = 0; j < nDim; j++ ) {
				if ( fabs(cell_momentum(icell,j)) > vel ) {
					vel = fabs(cell_momentum(icell,j));
				}
			}
			vel = vel/cell_gas_density(icell) + as;

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
