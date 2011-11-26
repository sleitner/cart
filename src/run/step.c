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
#include "system.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "gravity_step.h"
#include "hydro_step.h"
#include "hydro_tracer_step.h"
#include "particle_step.h"
#include "rt_step.h"
#include "starformation_step.h"
#include "starformation_feedback_step.h"
#include "step.h"


const int is_running = 1;

extern double t_init;
extern double t_end;

extern int min_time_refinement_factor;
extern int max_time_refinement_factor;
extern int time_refinement_level;

DEFINE_LEVEL_ARRAY(double,dtl);
DEFINE_LEVEL_ARRAY(double,dtl_old);
DEFINE_LEVEL_ARRAY(double,tl_old);

DEFINE_LEVEL_ARRAY(int,time_refinement_factor);
DEFINE_LEVEL_ARRAY(int,time_refinement_factor_old);

/*
//  The variable dt_global is the chosen global time-step without any 
//  extra reduction, and the variable frac_dt is the fraction of that 
//  time-step the code is going to take. After the code has recovered
//  from the CFL condition violation, frac_dt = 1.
*/
double dt_global, frac_dt;
int cfl_violation = 0;

#ifdef COSMOLOGY
extern double auni_init;
extern double auni_end;
DEFINE_LEVEL_ARRAY(double,abox_old);
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

extern double max_time_inc;
extern double min_time_dec;
extern double max_dt;
#ifdef COSMOLOGY
extern double max_a_inc;
extern double max_da;
#endif /* COSMOLOGY */

extern int steps_before_increasing;

int step_of_last_increase = 0;
int current_step_level = -1;

double min_courant_velocity = 1.0e-6;


void config_read_file(const char *filename);
void config_append_units_to_file(const char *filename);
void config_print_to_file(const char *filename, int append);
void config_plugins();

void init_run();
void run_output();

int global_timestep();
int timestep( int level, MPI_Comm level_com );

#ifdef HYDRO
void hydro_cfl_condition( int level, int *courant_cell, double *velocity );
#endif


void run( int restart, const char *restart_label ) {
	int current_steps;
	int level;

	init_logging( restart );

	if ( !restart ) {
		/* set up an initial decomposition */
		init_parallel_grid();
		init_tree();

		/* set up individual problem (responsible for setting
			time variables tl, dtl) on min_level only
		*/
		init_run();
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

#ifdef HYDRO 
	hydro_step_init();
#endif

	end_time( INIT_TIMER );

	current_steps = 0;
	last_restart_step = 0;

	dt_global = max_dt;
	frac_dt = 1.0;

#ifdef USER_PLUGIN
	config_plugins();
	PLUGIN_POINT(RunBegin)();
#else
	start_time( OUTPUT_TIMER );
	run_output();
	end_time( OUTPUT_TIMER );
#endif

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
				write_restart( WRITE_GENERIC, WRITE_GENERIC, WRITE_GENERIC, NULL );
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

		} else {
			step++;
			current_steps++;

#ifndef USER_PLUGIN
			if ( output_frequency > 0 && step % output_frequency == 0 ) {
				start_time( OUTPUT_TIMER );
				run_output();
				end_time( OUTPUT_TIMER );
			}
#endif
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

#ifdef USER_PLUGIN
	PLUGIN_POINT(RunEnd)();
#endif

	/* destroy buffer */
	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}
	
	cart_debug("total time = %f", end_time(TOTAL_TIME) );

	finalize_logging();
}


int global_timestep() {
	int level;
	int ret, global_ret;
	double fdt;

        cart_assert( buffer_enabled );

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
	ret = timestep( min_level, mpi.comm.run );
	current_step_level = -1;

	/* check if any other processors had problems (violation of CFL condition) */
	start_time( COMMUNICATION_TIMER );
	MPI_Allreduce( &ret, &global_ret, 1, MPI_INT, MPI_MIN, mpi.comm.run );

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
			modify( level, OP_REFINE );
		}
#endif /* REFINEMENT */

#ifdef STARFORM
#ifdef AGN               
                /* do agn mergers */
                agn_find_mergers();

        /* do agn formation */
#endif /* AGN */

		/* now remap ids of stars created in this timestep */
		remap_star_ids();
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER	
		rtStepEnd();
#endif /* RADIATIVE_TRANSFER */

#ifdef USER_PLUGIN
		PLUGIN_POINT(GlobalStepEnd)();
#endif

	} else {
		start_time( COMMUNICATION_TIMER );
		fdt = frac_dt;
		MPI_Allreduce( &fdt, &frac_dt, 1, MPI_DOUBLE, MPI_MIN, mpi.comm.run );
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
	int i, j, courant_cell;
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
	//  Backup hydro variables for updating fluxes on lower levels
	*/
	hydro_copy_vars( level, HYDRO_COPY_ALL );
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

	    cart_debug("Creating level %d...",level);
#ifdef COSMOLOGY
	    cart_debug("Spatial resolution: %le kpc (proper), %le CHIMP (comoving)",units->length*cell_size[level]/constants->kpc,units->length_in_chimps*cell_size[level]);
#else  /* COSMOLOGY */
	    cart_debug("Spatial resolution: %le kpc",units->length*cell_size[level]/constants->kpc);
#endif /* COSMOLOGY */


	    tl[level] = tl[level-1];
#ifdef COSMOLOGY
	    abox[level] = abox[level-1];
	    auni[level] = auni[level-1];
	    abox_old[level] = abox_old[level-1];
#endif /* COSMOLOGY */

	    num_steps_on_level[level] = num_steps_on_level[level-1]*time_refinement_factor[level];
	  }
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
	dt_needed = cfl_max * cell_size[level] / velocity;

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
                                max( cell_momentum(courant_cell,0), 
                                        max( cell_momentum(courant_cell,1), cell_momentum(courant_cell,2) )) /
                                cell_gas_density(courant_cell) ) * units->time / constants->Myr );
                cart_debug("---------------------------------------------------------");

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
		MPI_Allreduce( &ret, &true_ret, 1, MPI_INT, MPI_MIN, level_com);
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
#if defined(HYDRO) || defined(REFINEMENT)
    /* if hydro or refinement are enabled, we destroyed the value of the 
     * acceleration variable and need to recompute it here */
        compute_accelerations_particles(level);
#endif /* HYDRO || REFINEMENT */

        accelerate_particles(level);
#endif /* GRAVITY */

	move_particles( level );
	update_particle_list( level );

#ifdef STARFORM
        star_particle_feedback(level);

	/* update cell values changed by starformation and feedback */
	start_time( STELLAR_FEEDBACK_UPDATE_TIMER );
	update_buffer_level( level, all_hydro_vars, num_hydro_vars );
#ifdef AGN
	for ( i = level+1; i <= max_level; i++ ) {
                update_buffer_level( i, all_hydro_vars, num_hydro_vars );
        }
#endif /* AGN */
	end_time( STELLAR_FEEDBACK_UPDATE_TIMER );
#endif /* STARFORM */

        move_particles( level );
        update_particle_list( level );
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
	if(level < max_level)
	  {
	    modify( level, OP_REFINE | OP_DEREFINE );
	  }
#endif /* REFINEMENT */

#ifdef USER_PLUGIN
	PLUGIN_POINT(LevelStepEnd)(level,level_com);
#endif

        cart_debug("timestep(%u, %e %s, %d)", level, 
#ifdef COSMOLOGY
                dtl[level]*units->time/constants->Myr, "Myr",
#else
                dtl[level]*units->time/constants->s, "s",
#endif
                num_steps_on_level[level] );

	num_steps_on_level[level]++;

	end_time( LEVEL_TIMER );
	end_timing_level( level );

	return ret;
}


void set_timestepping_scheme()
{
  int level, i, j, lowest_level, synch_level;
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
  rtModifyTimeStep(&dt_local);
#endif /* RADIATIVE_TRANSFER */

  start_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
  MPI_Allreduce(&dt_local,&dt_new,1,MPI_DOUBLE,MPI_MIN,mpi.comm.run);
  end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );

#else 

  start_time( WORK_TIMER );

#ifdef HYDRO 
  hydro_cfl_condition(min_level,&courant_cell,&velocity);

  cart_debug("cfl cell %d [level %d]: velocity = %e cm/s, cs = %e cm/s, n = %e #/cc, dt = %e Myr", 
             courant_cell, cell_level(courant_cell), 
             max( cell_momentum(courant_cell,0),
                  max( cell_momentum(courant_cell,1),
                       cell_momentum(courant_cell,2) )) /
             cell_gas_density(courant_cell) * units->velocity/constants->cms,
             sqrt( cell_gas_gamma(courant_cell) * cell_gas_pressure(courant_cell) / 
                   cell_gas_density(courant_cell))*units->velocity/constants->cms,
             cell_gas_density(courant_cell)*units->number_density/constants->cc,
             cfl_run*cell_size[min_level]/velocity*units->time/constants->Myr );

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
  MPI_Allreduce(&dt_local,&dt_new,1,MPI_DOUBLE,MPI_MIN,mpi.comm.run);
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
  if(local_proc_id == MASTER_NODE) cart_debug("Ideal time-step: %lg s, previous: %lg s",
                dt_new*units->time/constants->s, dt_global*units->time/constants->s );
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
      if(local_proc_id == MASTER_NODE) cart_debug("Limited time-step increase to: %lg",dt_new*units->time/constants->s);
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
      if(local_proc_id == MASTER_NODE) cart_debug("Forced time-step decrease to: %lg",dt_new*units->time/constants->s);
#endif
    }

  end_time( WORK_TIMER );

#endif /* CONSTANT_TIMESTEP */

  start_time( WORK_TIMER );

  dtl[min_level] = dt_new;

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
  lowest_level = max_level_now_global(mpi.comm.run);

#ifdef HYDRO 
  for(level=min_level+1; level<=lowest_level; level++)
    {
      hydro_cfl_condition(level,&courant_cell,&velocity);

      if(velocity > 0.0)
	{
	  dt_new = cfl_run*cell_size[level]/velocity;
	  if(dt_new < dtl_local[level])
	    {
	      dtl_local[level] = dt_new;
	      cart_debug("new cfl cell %d [level %d]: velocity = %e cm/s, cs = %e cm/s, n = %e #/cc, dt = %e Myr", 
			 courant_cell, cell_level(courant_cell), 
			 max( cell_momentum(courant_cell,0), 
			      max( cell_momentum(courant_cell,1), cell_momentum(courant_cell,2) )) /
			 cell_gas_density(courant_cell) * units->velocity/constants->cms,
			 sqrt( cell_gas_gamma(courant_cell) * cell_gas_pressure(courant_cell) / 
			       cell_gas_density(courant_cell))*units->velocity/constants->cms,
			 cell_gas_density(courant_cell)*units->number_density/constants->cc,
			 cfl_run*cell_size[min_level]/velocity*units->time/constants->Myr );
	    }
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
      MPI_Allreduce(&dtl_local[min_level+1],&dtl[min_level+1],lowest_level-min_level,MPI_DOUBLE,MPI_MIN,mpi.comm.run);
      end_time( CHOOSE_TIMESTEP_COMMUNICATION_TIMER );
    }

  start_time( WORK_TIMER );

  /*
  // Synch time-steps on all levels above time_refinement_level or lowest_level.
  */
  synch_level = min(time_refinement_level,lowest_level);
  time_refinement_factor[synch_level] = 1;
  for(level=min_level; level<synch_level; level++)
    {
      time_refinement_factor[level] = 1;
      dtl[level] = dtl[synch_level];
    }

  /*
  // Satisfy min_time_refinement_factor limit on all levels below synch_level.
  */
  for(level=synch_level+1; level<=lowest_level; level++)
    {
      time_refinement_factor[level] = max(min_time_refinement_factor,(int)(0.999+dtl[level-1]/dtl[level]));
      dtl[level] = dtl[level-1]/time_refinement_factor[level];
    }

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
  //  their refinement scheme to the old-style, factor-of-two refinement
  //  (if it is allowed).
  */
  for(level=lowest_level+1; level<=max_level; level++)
    {
      time_refinement_factor[level] = min(2,max_time_refinement_factor);
      dtl[level] = dtl[level-1]/time_refinement_factor[level];
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
      cart_debug("chose %le %s as our next global time-step (%lg of the optimal value)",
#ifdef COSMOLOGY
		 dtl[min_level]*units->time/constants->Myr, "Myr",
#else
		 dtl[min_level]*units->time/constants->s, "s",
#endif
		 frac_dt );

      work = 1.0;
      for(level=min_level+1; level<=lowest_level; level++)
	{
	  work *= time_refinement_factor[level];
          cart_debug("level %d, dt = %9.3le %s, time-ref = %d, global time-ref = %lg", level,
#ifdef COSMOLOGY
		     dtl[level]*units->time/constants->Myr, "Myr",
#else
		     dtl[level]*units->time/constants->s, "s",
#endif
		     time_refinement_factor[level],work);
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

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

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
	cart_free( level_cells );
	
	*courant_cell = ivas;
	*velocity = vas;
}
#endif /* HYDRO */
