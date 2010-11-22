#include "config.h"

#include <limits.h>
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
#include "plugin.h"
#include "refinement.h"
#include "rt_solver.h"
#include "starformation.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


double t_init = -1.0e38;
double t_end = 1.0e38;
DEFINE_LEVEL_ARRAY(double,dtl);
DEFINE_LEVEL_ARRAY(double,dtl_old);
DEFINE_LEVEL_ARRAY(double,tl);
DEFINE_LEVEL_ARRAY(double,tl_old);

DEFINE_LEVEL_ARRAY(int,time_refinement_factor);
DEFINE_LEVEL_ARRAY(int,time_refinement_factor_old);

/*
//  Time refinement scheme defaults to the old factor-of-two refinement.
*/
int min_time_refinement_factor = 2;
int max_time_refinement_factor = 2;

/*
//  The variable dt_global is the chosen global time-step without any 
//  extra reduction, and the variable frac_dt is the fraction of that 
//  time-step the code is going to take. After the code has recovered
//  from the CFL condition violation, frac_dt = 1.
*/
double dt_global, frac_dt;
int cfl_violation = 0;

#ifdef COSMOLOGY
double auni_init = 1.0e-3;
double auni_end = 1.0;
DEFINE_LEVEL_ARRAY(double,abox);
DEFINE_LEVEL_ARRAY(double,abox_old);
DEFINE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

DEFINE_LEVEL_ARRAY(unsigned int,num_steps_on_level);

int max_steps = 0;
double timelimit = 0.0;

int max_mpi_sync_level = min_level;

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


int timestep( int level, MPI_Comm level_com );


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

  control_parameter_add3(control_parameter_int,&max_mpi_sync_level,"max-mpi-sync-level","max_mpi_sync_level","max_cfl_sync_level","maximum level at which MPI calls are synchronized across separate nodes.");

#ifdef HYDRO 
  control_parameter_add3(control_parameter_cflrun,&cfl_run,"cfl-run","cfl_run","cfl","the CFL number for setting the time-step.");

  control_parameter_add3(control_parameter_cflmax,&cfl_max,"cfl-max","cfl_max","cfl","the maximum acceptable CFL number. If this number is exceeded, the time-step needs to be redone. It is a good sense to set this number just a little bit higher than <cfl-run>, to avoid extra restarts due to numerical noise.");
#endif /* HYDRO */

#ifdef PARTICLES
  control_parameter_add2(control_parameter_double,&particle_cfl,"particle-cfl","particle_cfl","the CFL number for particle dynamics. In HYDRO mode this number is usually not needed, as the grid CFL conditions superceeds that of particles. Setting it to zero disables this limit.");
#endif /* PARTICLES */

  control_parameter_add2(control_parameter_double,&max_time_inc,"max-timestep-increment","max_time_inc","the largest factor by which the time-step is allowed to increase.");

  control_parameter_add2(control_parameter_double,&min_time_dec,"min-timestep-decrement","min_time_dec","the smallest factor by which the time-step is allowed to decrease.");

  control_parameter_add2(control_parameter_double,&max_dt,"max-dt","max_dt","maximum allowed time-step.");

  control_parameter_add2(control_parameter_int,&steps_before_increasing,"timesteps-before-increasing","steps_before_increasing","number of global time-steps to make before the time-step is allowed to increase.");

  control_parameter_add2(control_parameter_int,&min_time_refinement_factor,"time-refinement-factor:min","min_time_refinement_factor","minimum allowed refinement factor between the time-steps on two successive levels.");
  control_parameter_add2(control_parameter_int,&max_time_refinement_factor,"time-refinement-factor:max","max_time_refinement_factor","maximum allowed refinement factor between the time-steps on two successive levels. If <time-refinement-factor:min> = <time-refinement-factor:max>, then time refinement is done by a constant factor; <time-refinement-factor:min> = <time-refinement-factor:max> = 2 recovers the old factor-of-two scheme.");
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

  cart_assert(max_mpi_sync_level >= min_level);

#ifdef HYDRO 
  cart_assert(cfl_run>0.0 && !(cfl_run>cfl_max));
#endif /* HYDRO */

#ifdef PARTICLES
  cart_assert(!(particle_cfl < 0.0));
#endif /* PARTICLES */

  cart_assert(max_time_inc > 1.0);

  cart_assert(min_time_dec > 1.0);

  cart_assert(max_dt >= 0.0);

  cart_assert(steps_before_increasing > 0);

  cart_assert(min_time_refinement_factor > 0);
  cart_assert(max_time_refinement_factor >= min_time_refinement_factor);

#ifndef HYDRO
  if(min_time_refinement_factor==1 && particle_cfl==0.0)
    {
      cart_error("In a pure N-body mode, with <particle-cfl> parameter set to 0.0, it is incorrect to set the <time-refinement-factor:min> parameter to 1. In that case there is no condition to make time-steps on lower levels smaller than on the top level, so the particle trajectories will be integrated incorrectly.");
    }
#endif /* HYDRO */

  dt_global = max_dt;
  frac_dt = 1.0;
}


int global_timestep() {
	int level;
	int ret, global_ret;
	double fdt;

	start_time( LEVEL_TIMER );
	start_time( WORK_TIMER );

	current_step_level = -1;
	cfl_violation = 0;

	/* set old vars */
	for ( level = min_level; level <= max_level; level++ ) {
		tl_old[level] = tl[level];
		dtl_old[level] = dtl[level];
		time_refinement_factor_old[level] = time_refinement_factor[level];
#ifdef COSMOLOGY
		abox_old[level] = abox[level];
#endif /* COSMOLOGY */
		num_steps_on_level[level] = 0;
	}

	units_update(min_level);

	end_time( WORK_TIMER );

	/*
	//  Specify the hierarchy of time-steps
	*/
	set_timestepping_scheme();

#ifdef RADIATIVE_TRANSFER
	rtStepBegin();
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
	num_new_stars = 0;
	last_star_id = particle_species_indices[num_particle_species]-1;
#endif /* STARFORM */

#ifdef USER_PLUGIN
	PLUGIN_POINT(GlobalStepBegin)();
#endif
	ret = timestep( min_level, MPI_COMM_WORLD );
	current_step_level = -1;

#ifdef USER_PLUGIN
	PLUGIN_POINT(GlobalStepEnd)();
#endif

	/* check if any other processors had problems (violation of CFL condition) */
	start_time( COMMUNICATION_TIMER );
	MPI_Allreduce( &ret, &global_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );

	/*
	// NG: this should never happen!!! ret is already communicated to all interested parties!
	*/
	cart_assert(ret == global_ret);

	end_time( COMMUNICATION_TIMER );

	/* do a last refinement step (without allowing derefinement */
	if ( global_ret != -1 ) {
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		for ( level = min_level; level <= max_level-1; level++ ) {
			assign_density(level);
		}
#endif

#ifdef REFINEMENT
		for ( level = min_level; level <= max_level-1; level++ ) {
			modify( level, 0 );
		}
#endif

#ifdef STARFORM
		/* now remap ids of stars created in this timestep */
		remap_star_ids();
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER	
		rtStepEnd();
#endif /* RADIATIVE_TRANSFER */
	} else {
		start_time( COMMUNICATION_TIMER );
		fdt = frac_dt;
		MPI_Allreduce( &fdt, &frac_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
		end_time( COMMUNICATION_TIMER );
#ifdef DEBUG_TIMESTEP
		cart_debug("Reducing frac_dt to: %lg",frac_dt);
#endif
	}

	end_time(LEVEL_TIMER);
	return global_ret;
}

int timestep( int level, MPI_Comm level_com ) 
/* returns -1 if timestep would invalidate cfl condition */
{
	int j, courant_cell;
	double velocity;
	int ret;
	int step_ret;
	int true_ret;
	double dt_needed;
	MPI_Comm child_com;
	int refined;

	cart_assert( level >= min_level && level <= max_level );

	current_step_level = level;

	/* assume step was sucessful */
	ret = 0;

	start_timing_level( level );
	start_time( LEVEL_TIMER );
	
#ifdef USER_PLUGIN
	PLUGIN_POINT(LevelStepBegin)(level,level_com);
#endif

	start_time( LOWER_LEVEL_TIMER );  /* this is for internal accounting only */

	refined = (level<max_level && level<max_level_now());

        if(level <= max_mpi_sync_level)
          {
	    /* 
	    //  Create a child communicator
	    */
	    start_time( COMMUNICATION_TIMER );
	    MPI_Comm_split(level_com,refined,local_proc_id,&child_com);
	    end_time( COMMUNICATION_TIMER );
	  }
	else
	  {
	    child_com = level_com;
	  }

	start_time( WORK_TIMER );
	units_update(level);
	end_time( WORK_TIMER );


	/*
	//  NG: These two calls look like a waste - HOWEVER, they are needed
	//  for lower levels when they pick up values from the parent cells.
	//  Since we cannot restrict these calls to only refined nodes (they
	//  call update_buffer_level), we need to have them for all
	//  nodes, even if it is a waste for unrefined nodes.
	//
	//  Perhaps, in the future these calls may be restricted to only cells
	//  that border refined regions.
	*/
#ifdef HYDRO
	/*
	//  Backup hydro variables for upodating fluxes on lower levels
	*/
	hydro_copy_vars( level, COPY_ZERO_REF, COPY_SPLIT_NEIGHBORS );
#endif /* HYDRO */
#if defined(GRAVITY) && defined(PARTICLES)
	/*
	//  The accelerations on this level may be used by boundary
	//  cells on a lower level, so we need to pre-compute them.
	*/
	compute_accelerations_particles(level);
#endif /* GRAVITY && PARTICLES */

	if(refined)
	  {
	    for(j=0; j<time_refinement_factor[level+1]; j++)
	      {
	        step_ret = timestep(level+1,child_com);
		current_step_level = level;
		ret = min(ret,step_ret);
		if(ret==-1 && level<max_mpi_sync_level)
		  {
		    break;
		  }
	      }
	  }

        if(level <= max_mpi_sync_level)
          {
	    start_time( COMMUNICATION_TIMER );
	    MPI_Comm_free( &child_com );
	    end_time( COMMUNICATION_TIMER );
	  }

	end_time( LOWER_LEVEL_TIMER );

	/*
	//  Units need to be set here, since lower levels could have 
	//  reset the basic units to a further moment in time.
	*/
	start_time( WORK_TIMER );

	if(level>min_level && tl[level]+0.1*dtl[level]<tl[level-1])
	  {
	    /*
	    //  A new level appeared locally. Sync the time variables 
	    //  with the parent. The time-step variable should've been 
	    //  set already.
	    */
	    cart_assert(!refined);

	    cart_debug("Creating level %d",level);

	    tl[level] = tl[level-1];
#ifdef COSMOLOGY
	    abox[level] = abox[level-1];
	    auni[level] = auni[level-1];
	    abox_old[level] = abox_old[level-1];
#endif

	    num_steps_on_level[level] = num_steps_on_level[level-1]*time_refinement_factor[level];
	  }
	units_update(level);
	end_time( WORK_TIMER );

#ifdef HYDRO
	/* test if timestep is still valid */
	start_time( WORK_TIMER );
	hydro_cfl_condition( level, &courant_cell, &velocity );
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

		/*
		//  Set-up a time-step restriction: do it once only,
		//  otherwise consequitive arrivals to this place will
		//  reduce frac_dt to an unreasonably small value.
		*/
		if(!cfl_violation)
		  {
		    cfl_violation = 1;
		    frac_dt *= dt_needed/(min_time_dec*dtl[level]);
		  }
		ret = -1;
	}
	end_time( WORK_TIMER );

	if ( level <= max_mpi_sync_level ) {
		start_time( COMMUNICATION_TIMER );
		MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
		end_time( COMMUNICATION_TIMER );

		if ( true_ret < 0 ) {
#ifdef USER_PLUGIN
		        PLUGIN_POINT(LevelStepFail)(level,level_com);
#endif
			end_time( LEVEL_TIMER );
			end_timing_level( level ); 
			return true_ret;
		}
	}
#endif /* HYDRO */

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

#ifdef RADIATIVE_TRANSFER
	/* Do RT step */
	rtLevelUpdate(level);
	if(level <= max_mpi_sync_level)
	  {
	    rtGlobalUpdate(level,level_com);
	  }
#endif /* RADIATIVE_TRANSFER */

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
	compute_accelerations_particles(level);
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

        start_time( WORK_TIMER );
        units_update(level);
        end_time( WORK_TIMER );

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

#ifdef USER_PLUGIN
	PLUGIN_POINT(LevelStepEnd)(level,level_com);
#endif

	cart_debug("timestep(%u, %e, %d)", level, dtl[level], num_steps_on_level[level] );

	num_steps_on_level[level]++;

	end_time( LEVEL_TIMER );
	end_timing_level( level );

	return ret;
}


void set_timestepping_scheme()
{
  int level, i, j, lowest_level;
  int courant_cell, done;
  double velocity;
  double dt_local, dt_new, dda, work;
  DEFINE_LEVEL_ARRAY(double,dtl_local);

  start_time( CHOOSE_TIMESTEP_TIMER );

  /*
  //  Compute the top level step
  */
  dt_local = max_dt;

#ifdef CONSTANT_TIMESTEP

#ifdef RADIATIVE_TRANSFER
  rtModifyTimeStep(dt_local);
#endif /* RADIATIVE_TRANSFER */

  start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
  MPI_Allreduce(&dt_local,&dt_new,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );

#else 

  start_time( WORK_TIMER );

#ifdef HYDRO 
  hydro_cfl_condition(min_level,&courant_cell,&velocity);

  cart_debug("cfl cell: velocity = %e, cell = %u, pressure = %e, density = %e, momentum = %e %e %e",
	     velocity, courant_cell, cell_gas_pressure(courant_cell),
	     cell_gas_density(courant_cell), cell_momentum(courant_cell,0),
	     cell_momentum(courant_cell,1), cell_momentum(courant_cell,2) );

  if(velocity > 0.0)
    {
      dt_local = min(dt_local,cfl_run*cell_size[min_level]/velocity);
    }
#endif /* HYDRO */

#ifdef PARTICLES
  if(particle_cfl > 0.0)
    {
      velocity = 0.0;
      for(i=0; i<num_particles; i++) if(particle_level[i] == min_level)
	{
	  for(j=0; j<nDim; j++) velocity = max(fabs(particle_v[i][j]),velocity);
	}

      if(velocity > 0.0)
	{
	  dt_local = min(dt_local,particle_cfl*cell_size[min_level]/velocity);
	}
    }
#endif /* PARTICLES */

  end_time( WORK_TIMER );

  start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
  MPI_Allreduce(&dt_local,&dt_new,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );

  start_time( WORK_TIMER );

#ifdef COSMOLOGY
  /*
  // dt_new could be very large, it can break cosmology module; hence, combine with max_a_inc
  */
  dda = abox_from_tcode(min(tl[min_level]+dt_new,tcode_from_abox(abox[min_level]*max_a_inc))) - abox[min_level];

  if(max_da > 0.0) dda = min(dda,max_da);

  dt_new = min(dt_new,tcode_from_abox(abox[min_level]+dda)-tl[min_level]);

#endif /* COSMOLOGY */

  /*
  //  dt_new is the ideal time-step we should have.
  */

#ifdef DEBUG_TIMESTEP
  if(local_proc_id == MASTER_NODE) cart_debug("Ideal time-step: %lg, previous: %lg",dt_new,dt_global);
#endif

  /* 
  //  Enforce maximum change in timestep and frac_dt
  */
  if(dt_new > 1.001*frac_dt*dt_global) // allow for round-off error
    {
      if(step-step_of_last_increase >= steps_before_increasing)
	{
	  /*
	  //  Increase by at most max_time_inc
	  */
	  frac_dt *= max_time_inc;
	  if(frac_dt > 1.0)
	    {
	      dt_new = min(frac_dt*dt_global,dt_new);
	      frac_dt = 1.0;
	    }
	  step_of_last_increase = step;
	} 
      else
	{
	  dt_new = dt_global;
	}
#ifdef DEBUG_TIMESTEP
      if(local_proc_id == MASTER_NODE) cart_debug("Limited time-step increase to: %lg",dt_new);
#endif
    }
  else if(dt_new < 0.999*frac_dt*dt_global) // allow for round-off error
    {
      /*
      //  Decrease by at least min_time_dec (so that next time
      //  no decrease may be needed).
      */
      dt_new = min(dt_global/min_time_dec,dt_new);
      frac_dt = 1.0;
      step_of_last_increase = step; // forbid an increase right away
#ifdef DEBUG_TIMESTEP
      if(local_proc_id == MASTER_NODE) cart_debug("Forced time-step decrease to: %lg",dt_new);
#endif
    }

#endif /* CONSTANT_TIMESTEP */

  dtl[min_level] = dt_new;
  time_refinement_factor[min_level] = 1;

  /*
  //  Ensure consistency of time variables 
  */
  for(level=min_level+1; level<=max_level; level++)
    {
      tl[level] = tl[min_level];
      dtl_local[level] = dt_new;

#ifdef COSMOLOGY
      abox[level] = abox_from_tcode(tl[level]);
      auni[level] = auni_from_tcode(tl[level]);
#endif
    }

  /*
  //  Set the time refinements
  */
  lowest_level = max_level_now_global(MPI_COMM_WORLD);

#ifdef HYDRO 
  for(level=min_level+1; level<=lowest_level; level++)
    {
      hydro_cfl_condition(level,&courant_cell,&velocity);

      if(velocity > 0.0)
	{
	  dtl_local[level] = min(dtl_local[level],cfl_run*cell_size[level]/velocity);
	}
    }
#endif /* HYDRO */

#ifdef PARTICLES
  if(particle_cfl > 0.0)
    {
      for(i=0; i<num_particles; i++) if(particle_level[i]>min_level && particle_level[i]<=max_level)
	{
	  velocity = 0.0;
	  for(j=0; j<nDim; j++)	velocity = max(fabs(particle_v[i][j]),velocity);

	  if(velocity > 0.0)
	    {
	      dtl_local[particle_level[i]] = min(dtl_local[particle_level[i]],particle_cfl*cell_size[particle_level[i]]/velocity);
	    }
	}
    }
#endif /* PARTICLES */

  end_time( WORK_TIMER );

  if(lowest_level > min_level)
    {
      start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
      MPI_Allreduce(&dtl_local[min_level+1],&dtl[min_level+1],lowest_level-min_level,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
      end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
    }

  start_time( WORK_TIMER );

  for(level=min_level+1; level<=lowest_level; level++)
    {
      time_refinement_factor[level] = max(min_time_refinement_factor,1+(int)(dtl[level-1]/dtl[level]));
      dtl[level] = dtl[level-1]/time_refinement_factor[level];
    }

  /*
  // The time-stepping scheme is now done, but we are not guaranteed to satisfy max_time_refinement_factor limit.
  */
  do
    {
      done = 1;
      for(level=min_level+1; level<=lowest_level; level++)
	{
	  if(time_refinement_factor[level] > max_time_refinement_factor)
	    {
	      dda = (double)max_time_refinement_factor/time_refinement_factor[level];
	      for(j=min_level; j<level; j++)
		{
		  dtl[j] *= dda;
		}
	      time_refinement_factor[level] = max_time_refinement_factor;
	      done = 0;
	    }
	}
    }
  while(!done);

  /*
  //  Now we need to take care of levels that do not exist yet. Some of
  //  these levels can appear during the time-step. If their time-steps are
  //  not refined, the CFL condition is guaranteed to be violated as soon
  //  as a new level steps for the first time. So, as a first guess, set
  //  their refinement scheme to the old-style, factor-of-two refinement.
  */
  for(level=lowest_level+1; level<=max_level; level++)
    {
      time_refinement_factor[level] = 2;
      dtl[level] = 0.5*dtl[level-1];
    }

  /*
  //  Now we have a complete scheme without any extra reduction. Save it.
  */
  dt_global = dtl[min_level];

  /*
  //  Impose the reduction factor.
  */
  for(level=min_level; level<=max_level; level++)
    {
      dtl[level] *= frac_dt;
    }

  end_time( WORK_TIMER );

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("chose %le as our next global time-step (%lg of the optimal value)",dtl[min_level],frac_dt);

      work = 1.0;
      for(level=min_level+1; level<=lowest_level; level++)
	{
	  work *= time_refinement_factor[level];
	  cart_debug("level %d, dt = %9.3le, time-ref = %d, global time-ref = %lg",level,dtl[level],time_refinement_factor[level],work);
	}
    }

#ifdef STARFORM
  start_time( WORK_TIMER );
  /* compute frequency of star formation calls */
  work = 1.0;
  for(level=min_level; level<=max_level; level++)
    {
      work *= time_refinement_factor[level];
      star_formation_frequency[level] = max(1,nearest_int(min(work,sf_sampling_timescale*constants->yr/(units->time*dtl[level]))));
    }
  end_time( WORK_TIMER );
#endif /* STARFORM */

  end_time( CHOOSE_TIMESTEP_TIMER );
}


#ifdef HYDRO
void hydro_cfl_condition( int level, int *courant_cell, double *velocity ) {
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
