#include "config.h"

#include <math.h>
#include <string.h>

#include "agn.h"
#include "auxiliary.h"
#include "cell_buffer.h"
#include "cooling.h"
#include "cosmology.h"
#include "density.h"
#include "gravity.h"
#include "halo_finder.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "io.h"
#include "iterators.h"
#include "logging.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "plugin.h"
#include "refinement.h"
#include "rt.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "system.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "agn_step.h"
#include "cell_feedback_step.h"
#include "gravity_step.h"
#include "hydro_step.h"
#include "hydro_tracer_step.h"
#include "io_step.h"
#include "particle_step.h"
#include "rt_step.h"
#include "stardestruction_step.h"
#include "starformation_step.h"
#include "starformation_feedback_step.h"
#include "step.h"

extern int min_time_refinement_factor;
extern int max_time_refinement_factor;
extern int time_refinement_level;

#ifdef COSMOLOGY
extern double auni_end;
#else
extern double t_end;
#endif /* COSMOLOGY */

DEFINE_LEVEL_ARRAY(unsigned int,num_steps_on_level);

extern int max_steps;
extern double timelimit;

extern int max_mpi_sync_level;

#ifdef HYDRO 
extern double cfl_run;
extern double cfl_max;
#endif /* HYDRO */

#ifdef PARTICLES
extern double particle_cfl;
#endif /* PARTICLES */

extern double max_dt_inc;
extern double min_dt_dec;
extern double tol_dt_grow;
extern double max_dt;
#ifdef COSMOLOGY
extern double max_a_inc;
extern double max_da;
#endif /* COSMOLOGY */

/*
//  The factor to reduce the time-step with after a CFL violation. 
//  The time-step is reduced by a factor 1/(1+reduce_dt_factor), so
//  reduce_dt_factor=0 means no reduction.
*/
DEFINE_LEVEL_ARRAY(double,reduce_dt_factor);
DEFINE_LEVEL_ARRAY(double,dtl_try);
int post_cfl_violation = 0;

extern double reduce_dt_factor_shallow_dec;
extern double reduce_dt_factor_deep_dec;

extern int current_step_level;

double min_courant_velocity = 1.0e-6;

void config_read_file(const char *filename);
void config_append_units_to_file(const char *filename);
void config_print_to_file(const char *filename, int append);
void config_plugins();

void init_run();
void run_output();

int global_timestep();
int timestep( int level, MPI_Comm level_com );

void set_post_cfl_timestepping_scheme();
void satisfy_time_refinement_constraints(int lowest_level);
void finalize_timestepping_scheme();
#ifdef HYDRO
void hydro_cfl_condition( int level, int *courant_cell, double *velocity );
#endif


void run( int restart, const char *restart_label ) {
	int current_steps;
	int level;

	start_time( INIT_TIMER );

	init_logging( restart );

	if ( !restart ) {
		/* set up an initial decomposition */
		init_parallel_grid();
		init_tree();

		/* set up individual problem (responsible for setting
			time variables tl, dtl) on min_level only
		*/
		init_run();
		tl_old[min_level] = tl[min_level];
#ifdef COSMOLOGY
		abox[min_level] = abox_from_tcode(tl[min_level]);
		auni[min_level] = auni_from_tcode(tl[min_level]);
#endif /* COSMOLOGY */
		for ( level = min_level+1; level <= max_level; level++ ) {
			tl[level] = tl[min_level];
			tl_old[level] = tl[min_level];
#ifdef COSMOLOGY
			abox[level] = abox[min_level];
			auni[level] = auni[min_level];
#endif /* COSMOLOGY */
		}

		units_update(min_level);

		step = 0;
		current_output = 0;

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		for ( level = min_level; level <= max_level; level++ ) {
			cart_debug("assigning density on level %u", level );
			assign_density( level );

#ifdef GRAVITY
			cart_debug("solving potential on level %u", level );
			solve_poisson( level, 0 );
			cart_debug("done solving potential");

			if ( level > min_level+1 ) {
				cart_debug("restricting potential to level %u", level-1);
				restrict_to_level( level-1 );
			}
#endif		
		}
#endif /* GRAVITY || RADIATIVE_TRANSFER */

#ifdef HYDRO
		for ( level = min_level; level <= max_level; level++ ) {
			level_sweep_dir[level] = 0;
		}
#endif /* HYDRO */

#ifdef GRAVITY
#ifdef HYDRO
		for ( level = min_level; level <= max_level; level++ ) {
			copy_potential( level );
		}
#endif /* HYDRO */
#endif /* GRAVITY */
	} else {
		read_restart(restart_label);
		load_balance(); 

#ifdef COSMOLOGY
		abox[min_level] = abox_from_tcode(tl[min_level]);
		auni[min_level] = auni_from_tcode(tl[min_level]);
#endif /* COSMOLOGY */
		for ( level = min_level+1; level <= max_level; level++ ) {
			tl[level] = tl[min_level];
#ifdef COSMOLOGY
			abox[level] = abox[min_level];
			auni[level] = auni[min_level];
#endif /* COSMOLOGY */
		}

		for ( level = min_level; level <= max_level; level++ ) {
			cart_debug("num_cells_per_level[%u] = %u", level, num_cells_per_level[level] );
		}

#ifdef RADIATIVE_TRANSFER
		for ( level = min_level; level <= max_level; level++ ) {
			cart_debug("assigning density on level %u", level );
			assign_density( level );
		}
#endif /* RADIATIVE_TRANSFER */
	} 

#ifdef LOG_STAR_CREATION
        check_restart_star_creation();
        log_star_creation(-1,-1.0,FILE_OPEN);
#endif

#ifdef DEBUG
	check_map();
#endif

	config_append_units_to_file("config.log");
	config_print_to_file("history.log",1);

	end_time( INIT_TIMER );

	current_steps = 0;
	last_restart_step = 0;

	PLUGIN_POINT(RunBegin,());

	/*
	//  Plugin should use-up all un-extracted options
	*/
	die_on_unknown_options();

	for ( level = min_level; level <= max_level; level++ ) reduce_dt_factor[level] = 0.0;

	while ( 1 ) {

#ifdef COSMOLOGY
		cart_debug("a = %e, t = %e", auni[min_level], tl[min_level] );

		if ( auni[min_level] >= auni_end ) {
			break;
		}
#else
		cart_debug("t = %e", tl[min_level] );

		if ( tl[min_level] >= t_end ) {
			break;
		}
#endif /* COSMOLOGY */


		/*
		//  NG: I know it is a dangerous thing, but is very useful
		//  on large system, where the queue waiting times are long.
		*/
		if(file_exists("in-run-replace.cfg"))
		  {
		    config_read_file("in-run-replace.cfg");
		    config_print_to_file("config.log",0);
		  }

		if ( ( timelimit > 0.0 && current_steps > 0 && 
				timelimit-current_time(TOTAL_TIME,min_level-1) < 1.5*last_time(LEVEL_TIMER, min_level-1) ) ||
				( max_steps > 0 && current_steps >= max_steps ) ) {
			cart_debug("reached time or step limit... writing restart");

			/* avoid writing restart file twice on last step */
			if ( last_restart_step != step ) {
				start_time( RESTART_TIMER );
				destroy_cell_buffer();
				write_restart( WRITE_BACKUP, WRITE_BACKUP, WRITE_BACKUP );
				end_time( RESTART_TIMER );
			}

			/* do requeue command */
			if ( local_proc_id == MASTER_NODE && 
					strlen(requeue_command) > 0 ) {
				system( requeue_command );
			}
			break;
		}

		if ( global_timestep() == -1 ) {
			cart_debug("Error: could not complete timestep, restarting from previous timestep" );

			post_cfl_violation = 1;

			start_time( RESTART_TIMER );
			destroy_cell_buffer();
#ifdef PARTICLES
			init_particles();
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
			init_hydro_tracers();
#endif /* HYDRO_TRACERS */

			/*
			//  If we started from a labeled file and encountered
			//  a CFL violation in the first step, we need to 
			//  restart from the same label.
			*/
			read_restart(restart_label);
			load_balance(); 
			end_time( RESTART_TIMER );

#ifdef RADIATIVE_TRANSFER
			for ( level = min_level; level <= max_level; level++ ) {
			  cart_debug("assigning density on level %u", level );
			  assign_density( level );
			}
#endif /* RADIATIVE_TRANSFER */
	
#ifdef LOG_STAR_CREATION
			cart_debug("wipe temp_star files in case of CFL restart");
			log_star_creation(-1,-1.0,FILE_CLOSE); 
			log_star_creation(-1,-1.0,FILE_OPEN); 
#endif
			PLUGIN_POINT(CFLRestart,());
		} else {
			current_steps++;

			post_cfl_violation = 0;

			save_check();
			/*
			//  If we started from a labeled file and wrote the 
			//  the restart file, then we need to erase the 
			//  restart_label, so that we start from a written
			//  restart file.
			*/
			if(last_restart_step == step) restart_label = NULL;

			if ( load_balance_frequency > 0 && step % load_balance_frequency == 0 ) {
				load_balance();
			}

			cart_debug("done with timestep %u at tl = %e", step, tl[min_level] );
			start_time( IO_TIMER );
			log_diagnostics();
			end_time( IO_TIMER );
		}
	}

	PLUGIN_POINT(RunEnd,());

	/* destroy buffer */
	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}
	
	cart_debug("total time = %f", current_time(TOTAL_TIME,min_level-1) );

	finalize_logging();
}


int global_timestep() {
	int level;
	int ret, global_ret;
#if defined(COSMOLOGY) && defined(PARTICLES)
	halo_list *halos;
#endif

	cart_assert( buffer_enabled );

	start_time( LEVEL_TIMER );
	start_time( WORK_TIMER );

	current_step_level = -1;

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
	//  Specify the hierarchy of time-steps. We can come here in two cases: as a new step,
	//  or after va CFL violation. In the latter case we don't care about any time-stepping
	//  scheme, we should just reduce the last step further.
	*/
	if(post_cfl_violation == 0)
	  {
	    set_timestepping_scheme();
	  }
	else
	  {
	    set_post_cfl_timestepping_scheme();
	  }

	/*
	// Reset reduction factors for the new step, save the step we
	// are trying
	*/
	for(level=min_level; level<=max_level; level++)
	  {
	    reduce_dt_factor[level] = 0.0;
	    dtl_try[level] = dtl[level];
	  }

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

#ifdef STAR_FORMATION
	start_star_allocation();	
#endif /* STAR_FORMATION */

	PLUGIN_POINT(GlobalStepBegin,());

	ret = timestep( min_level, mpi.comm.run );
	current_step_level = -1;
	step++;

	/* check if any other processors had problems (violation of CFL condition) */
	start_time( COMMUNICATION_TIMER );
	MPI_Allreduce( &ret, &global_ret, 1, MPI_INT, MPI_MIN, mpi.comm.run );

	/*
	// NG: this should never happen!!! ret is already communicated to all interested parties!
	*/
	cart_assert(ret == global_ret);

	end_time( COMMUNICATION_TIMER );

	if ( global_ret != -1 ) {
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		for ( level = min_level; level <= max_level-1; level++ ) {
			assign_density(level);
		}
#endif

#ifdef REFINEMENT
		for ( level = min_level; level <= max_level-1; level++ ) {
			modify( level, OP_REFINE );
		}
#endif /* REFINEMENT */

#if defined(COSMOLOGY) && defined(PARTICLES)
		if ( halo_finder_frequency > 0 && step % halo_finder_frequency == 0 ) {
			halos = find_halos();
			if ( halos->num_halos > 0 ) {
				write_halo_list(halos);
			}
		} else {
			halos = NULL;
		}

#ifdef STAR_FORMATION
#ifdef AGN               
		/* do agn mergers */
		agn_find_mergers();

        /* do agn seeding */
		agn_seed( halos );
#endif /* AGN */

		/* now remap ids of stars created in this timestep */
		start_time( COMMUNICATION_TIMER );
		end_star_allocation();
		end_time( COMMUNICATION_TIMER );
#endif /* STAR_FORMATION */

		if ( halos != NULL ) {
			destroy_halo_list(halos);
		}
#endif /* COSMOLOGY && PARTICLES */

#ifdef RADIATIVE_TRANSFER	
		rtStepEnd();
#endif /* RADIATIVE_TRANSFER */

		PLUGIN_POINT(GlobalStepEnd,());

	}

	end_time(LEVEL_TIMER);
	return global_ret;
}

int timestep( int level, MPI_Comm level_com ) 
/* returns -1 if timestep would invalidate cfl condition */
{
	int j, courant_cell;
	double velocity;
	int nlevel;
	int factor;
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

	PLUGIN_POINT(LevelStepBegin,(level,level_com));

	start_time( WORK_TIMER );
	units_update(level);
	end_time( WORK_TIMER );

	/*
	//  NG: These two calls look like a waste - HOWEVER, they are needed
	//  for lower levels when they pick up values from the parent cells.
	//  Since we cannot restrict these calls to only refined nodes (they
	//  call update_buffer_level), we need to have them for all
	//  nodes, even if it is a waste for unrefined nodes.
	*/
#ifdef HYDRO
	/*
	//  Backup hydro variables for updating fluxes on lower levels and 
	//   for initial hydro level step
	*/
	hydro_copy_vars( level, HYDRO_COPY_ALL );
#endif /* HYDRO */
#if defined(GRAVITY) && defined(PARTICLES)
	/*
	//  The accelerations on this level may be used by boundary
	//  cells on a lower level, so we need to pre-compute them.
	//  Also necessary even if !refined, due to update_buffer_level
	//  call.
	*/
	compute_accelerations_particles(level);
#endif /* GRAVITY && PARTICLES */

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

	if(refined)
	{
		for(j=0; j<time_refinement_factor[level+1]; j++)
		{
			step_ret = timestep(level+1,child_com);
			current_step_level = level;
			ret = MIN(ret,step_ret);
			if(ret==-1 && level<max_mpi_sync_level)
			{
				break;
			}
		}
	} else {
		/* advance timestep on lower levels */
		factor = 1;
		for ( nlevel = level + 1; nlevel <= max_level; nlevel++ ) {
			tl[nlevel] = tl[level] + dtl[level];
			tl_old[nlevel] = tl[nlevel] - dtl[nlevel];

#ifdef COSMOLOGY
			abox[nlevel] = abox_from_tcode( tl[nlevel] );
			auni[nlevel] = auni_from_tcode( tl[nlevel] );
			abox_old[nlevel] = abox_from_tcode( tl_old[nlevel] );
#endif

			factor *= time_refinement_factor[nlevel];
			num_steps_on_level[nlevel] += factor;
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
	units_update(level);
	end_time( WORK_TIMER );

#ifdef HYDRO

#ifdef BLASTWAVE_FEEDBACK
	/* check time precision ... */
	check_bwtime_precision(level);
#endif /* BLASTWAVE_FEEDBACK */

	/* test if timestep is still valid */
	start_time( WORK_TIMER );
	hydro_cfl_condition( level, &courant_cell, &velocity );

	/* velocity can be 0 if this level has no leaves */
	if(velocity > 0.0)
	  {
	    dt_needed = cfl_max * cell_size[level] / velocity;
	  }
	else
	  {
	    dt_needed = dtl[level];
	  }

	/* check for cfl condition violation... */
	if ( dtl[level] > dt_needed && ret != -1 ) {
		cart_debug("---------------------------------------------------------");
		cart_debug("CFL CONDITION VIOLATED:");
		cart_debug("current dt = %.25e Myr", dtl[level]*units->time / constants->Myr );
		cart_debug("needed  dt = %.25e Myr", dt_needed*units->time / constants->Myr );
		cart_debug("CFL tolerance = %.3f / %.3f = %.4e", cfl_max, cfl_run, cfl_max/cfl_run );
		cart_debug("courant cell information:");
		cart_debug("T  = %e K", cell_gas_temperature(courant_cell)*units->temperature/constants->K );
		cart_debug("P  = %e ergs cm^-3", cell_gas_pressure(courant_cell)*units->energy_density/constants->barye );
		cart_debug("n  = %e cm^-3", cell_gas_density(courant_cell)*units->number_density*constants->cc );
		cart_debug("cs = %e cm/s", sqrt( cell_gas_gamma(courant_cell) * cell_gas_pressure(courant_cell) / 
					cell_gas_density(courant_cell))*units->velocity/constants->cms );
		cart_debug("v  = %e %e %e cm/s", 
				cell_momentum(courant_cell,0)/cell_gas_density(courant_cell)*units->velocity/constants->cms, 
				cell_momentum(courant_cell,1)/cell_gas_density(courant_cell)*units->velocity/constants->cms,
				cell_momentum(courant_cell,2)/cell_gas_density(courant_cell)*units->velocity/constants->cms );
		cart_debug("dt sound crossing = %e Myr", ( cell_size[level] / 
					sqrt( cell_gas_gamma(courant_cell) * cell_gas_pressure(courant_cell) / cell_gas_density(courant_cell)) )*
				units->time / constants->Myr );
		cart_debug("dt bulk velocity  = %e Myr", 
				cell_size[level]/(
					MAX( cell_momentum(courant_cell,0), 
						MAX( cell_momentum(courant_cell,1), cell_momentum(courant_cell,2) )) /
					cell_gas_density(courant_cell) ) * units->time / constants->Myr );
		cart_debug("---------------------------------------------------------");

		/*
		//  Set-up a time-step restriction (but only the first time,
		//  other times the solution is already bogus).
		*/
		if(reduce_dt_factor[level] == 0.0)
		{
			/*
			//  It does not make much sense for rdf to be too small, it may get swallowed 
			//  by the changes in time-steps due to time-refinement being done by integer factors
			*/
			reduce_dt_factor[level] = MAX(min_dt_dec,dtl[level]/MAX(0.1*dtl[level],dt_needed)) - 1.0;
		}
		ret = -1;
	}
	end_time( WORK_TIMER );

	if ( level <= max_mpi_sync_level ) {
		start_time( COMMUNICATION_TIMER );
		MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, level_com);
		end_time( COMMUNICATION_TIMER );

		if ( true_ret < 0 ) {
			PLUGIN_POINT(LevelStepFail,(level,level_com));
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
#ifdef STAR_FORMATION
	if ( num_steps_on_level[level] % star_formation_frequency[level] == 0 ) {
		star_formation( level, star_formation_frequency[level] );
	}
#endif /* STAR_FORMATION */
#endif /* HYDRO */

#ifdef GRAVITY
#if defined(HYDRO) || defined(REFINEMENT)
	/* if hydro or refinement are enabled, we destroyed the value of the 
     * acceleration variable and need to recompute it here */
	compute_accelerations_particles(level);
#endif /* HYDRO || REFINEMENT */

	accelerate_particles(level);
#endif /* GRAVITY */

#ifdef STAR_FORMATION
	if ( num_steps_on_level[level] % star_feedback_frequency[level] == 0 ) {
		star_particle_feedback(level, star_feedback_frequency[level]);

#if defined(STAR_PARTICLE_TYPES) && defined(AGN)
		agn_feedback( level, star_feedback_frequency[level] );
#endif /* STAR_PARTICLE_TYPES && AGN */

		cell_feedback(level, star_feedback_frequency[level]);

#ifdef STAR_PARTICLE_TYPES 
		/* feedback should care about star destruction if it happens before t_next */
		star_destruction( level ); 
#endif /* STAR_PARTICLE_TYPES */

		/* update cell values changed by starformation and feedback */
		start_time( STELLAR_FEEDBACK_UPDATE_TIMER );
		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
#ifdef AGN
		for ( j = level+1; j <= max_level; j++ ) {
			update_buffer_level( j, all_hydro_vars, num_hydro_vars );
		}
#endif /* AGN */
		end_time( STELLAR_FEEDBACK_UPDATE_TIMER );
	}  
#endif /* STAR_FORMATION */

	move_particles( level );
	update_particle_list( level );
#endif /* PARTICLES */

#ifdef ISOTROPIC_TURBULENCE_ENERGY 
	hydro_isotropic_turbulence_sources(level);
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

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
	if(level < max_level)
	  {
	    modify( level, OP_REFINE | OP_DEREFINE );
	  }
#endif /* REFINEMENT */

        cart_debug("timestep(%u, %9.3e %s, %d)", level, 
#ifdef COSMOLOGY
				dtl[level]*units->time/constants->Myr, "Myr",
#else
                dtl[level]*units->time, "s",
#endif
                num_steps_on_level[level] );

	PLUGIN_POINT(LevelStepEnd,(level,level_com));

	num_steps_on_level[level]++;

	end_time( LEVEL_TIMER );
	end_timing_level( level );

	return ret;
}


#ifdef CONSTANT_TIMESTEP

void set_timestepping_scheme()
{
  int level;
  double dt_new;

  start_time( CHOOSE_TIMESTEP_TIMER );
  start_time( WORK_TIMER );

  dt_new = max_dt;
#ifdef RADIATIVE_TRANSFER
  rtModifyTimeStep(&dt_new);
#endif /* RADIATIVE_TRANSFER */

  for(level=min_level; level<=max_level; level++)
    {
      tl[level] = tl[min_level];
      time_refinement_factor[level] = 2;
      dtl[level] = dt_new*pow(0.5,level-min_level);
    }

  end_time( WORK_TIMER );

  finalize_timestepping_scheme();

  end_time( CHOOSE_TIMESTEP_TIMER );
}

#else  /* CONSTANT_TIMESTEP */

void set_timestepping_scheme()
{
  int level, i, j, lowest_level;
  int courant_cell;
  double velocity;
  double dt_new, dda;
  DEFINE_LEVEL_ARRAY(double,dtl_local);

  start_time( CHOOSE_TIMESTEP_TIMER );
  start_time( WORK_TIMER );

  /*
  //  Ensure consistency of time variables (just in case)
  */
  for(level=min_level+1; level<=max_level; level++)
    {
      tl[level] = tl[min_level];

#ifdef COSMOLOGY
      abox[level] = abox_from_tcode(tl[level]);
      auni[level] = auni_from_tcode(tl[level]);
#endif
    }

  /*
  //  Now compute the best CFL-limited time-step on each level,
  //  ignoring any restriction for now (we add them later).
  //  Only use levels that actually exist.
  */
  lowest_level = max_level_now_global(mpi.comm.run);

  /*
  //  This is just in case we have no CFL restrictions - something must be
  //  set as a time-step on each level
  */
  for(level=min_level; level<=max_level; level++)
      {
          if(max_dt > 0)
              dtl_local[level] = max_dt;
          else
              dtl_local[level] = 1e30;
      }

#ifdef HYDRO 
  for(level=min_level; level<=lowest_level; level++)
    {
      hydro_cfl_condition(level,&courant_cell,&velocity);

      if(velocity > 0.0)
	{
	  dt_new = cfl_run*cell_size[level]/velocity;
	  dtl_local[level] = MIN(dtl_local[level],dt_new);
	  cart_debug("cfl cell %d [level %d]: velocity = %e cm/s, cs = %e cm/s, n = %e #/cc, dt = %e Myr", 
		     courant_cell, cell_level(courant_cell), 
		     velocity * units->velocity/constants->cms,
		     sqrt( cell_gas_gamma(courant_cell) * cell_gas_pressure(courant_cell) / 
			   cell_gas_density(courant_cell))*units->velocity/constants->cms,
		     cell_gas_density(courant_cell)*units->number_density/constants->cc,
		     cfl_run*cell_size[min_level]/velocity*units->time/constants->Myr );
	}
    }
#endif /* HYDRO */

#ifdef PARTICLES
  if(particle_cfl > 0.0)
    {
      for(i=0; i<num_particles; i++) if(particle_level[i]>=min_level && particle_level[i]<=max_level)
	{
	  velocity = 0.0;
	  for(j=0; j<nDim; j++)	velocity = MAX(fabs(particle_v[i][j]),velocity);

	  if(velocity > 0.0)
	    {
	      dtl_local[particle_level[i]] = MIN(dtl_local[particle_level[i]],particle_cfl*cell_size[particle_level[i]]/velocity);
	    }
	}
    }
#endif /* PARTICLES */

  end_time( WORK_TIMER );

  start_time( COMMUNICATION_TIMER );
  start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
  MPI_Allreduce(&dtl_local[min_level],&dtl[min_level],lowest_level-min_level+1,MPI_DOUBLE,MPI_MIN,mpi.comm.run);
  end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
  end_time( COMMUNICATION_TIMER );

  start_time( WORK_TIMER );

#ifdef COSMOLOGY
  /*
  // If we have cosmological constraints, satisfy them too.
  */
  dda = abox_from_tcode(MIN(tl[min_level]+dtl[min_level],tcode_from_abox(abox[min_level]*max_a_inc))) - abox[min_level];
  if(max_da > 0.0) dda = MIN(dda,max_da);

  dtl[min_level] = MIN(dtl[min_level],tcode_from_abox(abox[min_level]+dda)-tl[min_level]);
#endif /* COSMOLOGY */

  /*
  //  Do what the function says it does...
  */
  satisfy_time_refinement_constraints(lowest_level);
  dt_new = dtl[min_level];

  /*
  //  For the first step, we are done. For later steps, however, we need to ensure that
  //  the time-step does not grow too fast or does not decrease too litle.
  */
  if(step > 0) for(level=min_level; level<=lowest_level; level++)
    {
      cart_assert(dtl_old[level] > 0.0);

      /*
      //  It is not possible to maintain the old behavior - to allow increase
      //  only after a given number of steps, because different levels can increase
      //  and descrease their time-steps simultaneously. But an even better strategy
      //  is to allow increase ONLY if the previous step is way too small - in real
      //  simulations the time-step usually only decreases with time.
      //
      //  How much is too small? Say, twice below min_time_dec
      */
      if(dtl[level] > tol_dt_grow*dtl_old[level])
	{
	  dtl[level]  = MIN(dtl[level],dtl_old[level]*max_dt_inc);
#ifdef DEBUG_TIMESTEP
	  if(local_proc_id == MASTER_NODE) cart_debug("Limiting the time-step increase at level %d to: %lg Myr",level,dtl[level]*units->time/constants->Myr);
#endif
	}
      else if(dtl[level] < 0.999*dtl_old[level]) /* allow for round-off error */
	{
	  /*
	  //  Decrease by at least min_dt_dec (so that next time no decrease may be needed).
	  */
	  dtl[level] = MIN(dtl[level],dtl_old[level]/min_dt_dec);
#ifdef DEBUG_TIMESTEP
	  if(local_proc_id == MASTER_NODE) cart_debug("Extending the time-step decrease at level %d to: %lg Myr",level,dtl[level]*units->time/constants->Myr);
#endif
	}
      else
	{
	  /*
	  //  Try to keep the step constant
	  */
	  dtl[level] = dtl_old[level];
	}
    }

  /*
  //  These changes could have violated the constraints, so impose them again.
  */
  satisfy_time_refinement_constraints(lowest_level);

  /*
  //  Now we need to take care of levels that do not exist yet. Some of
  //  these levels can appear during the time-step. If their time-steps are
  //  not refined, the CFL condition is guaranteed to be violated as soon
  //  as a new level steps for the first time. So, as a first guess, set
  //  their refinement scheme to the old-style, factor-of-two refinement
  //  (if it is allowed).
  */
  for(level=lowest_level+1; level<=max_level; level++)
    {
      time_refinement_factor[level] = MIN(2,max_time_refinement_factor);
      dtl[level] = dtl[level-1]/time_refinement_factor[level];
    }

  end_time( WORK_TIMER );

  finalize_timestepping_scheme();

  end_time( CHOOSE_TIMESTEP_TIMER );
}

#endif /* CONSTANT_TIMESTEP */


void set_post_cfl_timestepping_scheme()
{
  int level, j;
  DEFINE_LEVEL_ARRAY(double,fdt);

  start_time( CHOOSE_TIMESTEP_TIMER );

  /*
  //  Find the global time-step reduction factors
  */
  start_time( COMMUNICATION_TIMER );

  for(level=min_level; level<=max_level; level++)
    {
      fdt[level] = reduce_dt_factor[level];
    }
		   
  MPI_Allreduce(&fdt[min_level],&reduce_dt_factor[min_level],max_level-min_level+1,MPI_DOUBLE,MPI_MAX,mpi.comm.run);

  end_time( COMMUNICATION_TIMER );

#ifdef DEBUG_TIMESTEP
  for(level=min_level; level<=max_level; level++)
    {
      if(reduce_dt_factor[level] > 0.0)
	{
	  cart_debug("Must reduce timestep[%d] by factor=%lg",level,reduce_dt_factor[level]+1);
	}
    }
#endif

  start_time( WORK_TIMER );

  /*
  //  Impose the reduction factor(s). We do that in a complex way,
  //  by reducing the offending level and levels around it, but the latter with progressively
  //  lower reduction factors. We extend the process into the non-existent levels as
  //  well, just in case. New reduction factors go into fdt array.
  */
  for(level=min_level; level<=max_level; level++)
    {
      fdt[level] = 0.0;
    }

  for(level=min_level; level<=max_level; level++)
    {
      if(reduce_dt_factor[level] > 0.0)
	{
	  fdt[level] = MAX(fdt[level],reduce_dt_factor[level]);  /* we may have been already reduced by some other level */
	  for(j=min_level; j<level; j++)
	    {
	      fdt[j] = MAX(fdt[j],reduce_dt_factor[level]*pow(reduce_dt_factor_shallow_dec,level-j));
	    }
	  for(j=level+1; j<=max_level; j++)
	    {
	      fdt[j] = MAX(fdt[j],reduce_dt_factor[level]*pow(reduce_dt_factor_deep_dec,j-level));
	    }
	}
    }

  /*
  //  dtl was overwritten on restart, but dtl_try contains the actual steps which were taken
  */
  for(level=min_level; level<=max_level; level++)
    {
#ifdef DEBUG_TIMESTEP
      if(fdt[level] > 0.0)
	{
	  cart_debug("reducing dt at level[%d]: %12.6le -> %12.6le, fdt=%lg",level,dtl[level],dtl[level]/(1+fdt[level]),fdt[level]);
	}
#endif
      dtl[level] = dtl_try[level]/(1+fdt[level]);
    }

  /*
  //  These changes could have violated the constraints, so impose them again.
  */
  satisfy_time_refinement_constraints(max_level);

  end_time( WORK_TIMER );

  finalize_timestepping_scheme();

  end_time( CHOOSE_TIMESTEP_TIMER );
}


void satisfy_time_refinement_constraints(int lowest_level)
{
  int level, synch_level, done, j;
  double q;

  /*
  // Synch time-steps on all levels above time_refinement_level or lowest_level.
  */
  synch_level = MIN(time_refinement_level,lowest_level);

  time_refinement_factor[synch_level] = 1;
  for(level=min_level; level<synch_level; level++)
    {
      time_refinement_factor[level] = 1;
      dtl[level] = dtl[synch_level];
    }

  /*
  // Satisfy min_time_refinement_factor limit on all levels below synch_level.
  */
#ifdef DEBUG_TIMESTEP
  DEFINE_LEVEL_ARRAY(double,dt0);
  for(level=min_level; level<=lowest_level; level++) dt0[level] = dtl[level];
#endif

  for(level=synch_level+1; level<=lowest_level; level++)
    {
      time_refinement_factor[level] = MAX(min_time_refinement_factor,(int)(0.999+dtl[level-1]/dtl[level]));
      dtl[level] = dtl[level-1]/time_refinement_factor[level];
    }

#ifdef DEBUG_TIMESTEP
  DEFINE_LEVEL_ARRAY(double,dt1);
  for(level=min_level; level<=lowest_level; level++) dt1[level] = dtl[level];
#endif

  /*
  // Satisfy max_time_refinement_factor limit on all levels below synch_level.
  */
  do
    {
      done = 1;
      for(level=synch_level+1; level<=lowest_level; level++)
	{
	  if(time_refinement_factor[level] > max_time_refinement_factor)
	    {
	      q = (double)max_time_refinement_factor/time_refinement_factor[level];
	      for(j=min_level; j<level; j++)
		{
		  dtl[j] *= q;
		}
	      time_refinement_factor[level] = max_time_refinement_factor;
	      cart_assert(fabs(dtl[level]/dtl[level-1]*time_refinement_factor[level]-1.0) < 1.0e-5);
	      done = 0;
	    }
	}
    }
  while(!done);

#ifdef DEBUG_TIMESTEP
  cart_debug("satisfy_time_refinement_constraints:");
  for(level=min_level; level<=lowest_level; level++)
    {
      cart_debug("level[%d]: %12.6le -> %12.6le -> %12.6le",level,dt0[level],dt1[level],dtl[level]);
    }
#endif
}


void finalize_timestepping_scheme()
{
  int level, lowest_level;
  double work;

  lowest_level = max_level_now_global(mpi.comm.run);

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("chose %le %s as our next global time-step",
#ifdef COSMOLOGY
		 dtl[min_level]*units->time/constants->Myr, "Myr"
#else
		 dtl[min_level]*units->time, "s"
#endif
		 );

      work = 1.0;
      for(level=min_level+1; level<=lowest_level; level++)
	{
	  work *= time_refinement_factor[level];
          cart_debug("level %d, dt = %9.3le %s, time-ref = %d, global time-ref = %lg", level,
#ifdef COSMOLOGY
		     dtl[level]*units->time/constants->Myr, "Myr",
#else
		     dtl[level]*units->time, "s",
#endif
		     time_refinement_factor[level],work);
	}
    }

#ifdef STAR_FORMATION
  start_time( WORK_TIMER );
  /* compute frequency of star formation calls & star feedback calls */
  work = 1.0;
  for(level=min_level; level<=max_level; level++)
    {
      work *= time_refinement_factor[level];
      star_formation_frequency[level] = MAX(1,nearest_int(MIN(work,sf_sampling_timescale*constants->yr/(units->time*dtl[level]))));
      star_feedback_frequency[level] = MAX(1,nearest_int(MIN(work,feedback_sampling_timescale*constants->yr/(units->time*dtl[level]))));
    }
  end_time( WORK_TIMER );
#endif /* STAR_FORMATION */
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

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

                as = cell_gas_sound_speed(icell) ;
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
	cart_free( level_cells );
	
	*courant_cell = ivas;
	*velocity = vas;
}
#endif /* HYDRO */
