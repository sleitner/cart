#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "defs.h"
#include "config.h"
#include "io.h"
#include "tree.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
#include "refinement.h"
#include "refinement_operations.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "gravity.h"
#include "density.h"
#include "logging.h"
#include "auxiliary.h"
#include "cooling.h"

#ifdef _OPENMP
#include <omp.h>
#include <sys/utsname.h>
#endif

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#endif


void init_run();
void run_output();

int num_options = 0;
char **options = NULL;


int main ( int argc, char *argv[]) {
	int i;
	int current_steps;
	int restart;
	int level;
	double dt, restart_dt;
	double restart_a;
	const char *tmp;
	char c;

	#ifdef _OPENMP
	struct utsname uname_info;
	#endif

	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &local_proc_id );
	
	if ( num_procs > MAX_PROCS ) {
		cart_error("Number of processors exceeds limit! (%u > %u) The executable must be recompiled.\n", 
			num_procs, MAX_PROCS );
	}

	/* load configuration file */
	if ( argc < 2 ) {
		cart_error("Usage: art <config_file> [ restart_flag [ options ] ]");
	}

	read_config( argv[1] );	

	if ( argc == 2 ) {
		restart = 0;
		/* skip config file name */
		num_options = argc - 2;
		options = argv + 2;
	} else {
		/*
		//  Also support an option-style restart in the form
		//    -restart[=<value>]
		//  where <value> is the value of the scale factor that
		//  labels restart files.
		*/
		tmp = check_option1(argv[2],"restart","last");
		if(tmp != NULL)
		  {
		    restart = 1;
		    if(strcmp(tmp,"last") == 0)
		      {
			restart_a = 0.0;
		      }
		    else
		      {
			/*
			//  Read in restart value
			*/
			if(sscanf(tmp,"%lg",&restart_a) != 1)
			  {
			    cart_error("Invalid option value %s",argv[2]+9);
			  }
		      }
		  }
		else
		  {
		    if(sscanf(argv[2],"%d%c",&restart,&c) != 1)
		      {
			cart_error("Invalid restart flag %s\nValid values are either an integer number or an option '-restart[=<value>]'",argv[2]);
		      }
		    restart_a = 0.0;
		  }
		/* skip config file name and restart flag */
		num_options = argc - 3;
		options = argv + 3;
	#ifdef _OPENMP
		if ( argc > 3 )
		  {
			tmp = check_option1(argv[3],"omp",NULL);
		    if(tmp != NULL)
		      {
			if(sscanf(tmp,"%d",&i)!=1 || i<1 || i>256)
			  {
			    cart_error("Invalid option value %s",argv[3]+4);
			  }
			omp_set_num_threads(i);
			num_options--;
			options++;
		      }
		  }
	#endif
	}

	cart_debug("my local pid = %u", getpid() );

	#ifdef _OPENMP
	cart_debug("num openmp threads = %u", omp_get_max_threads() );

	if ( num_procs > 1 ) {
		uname(&uname_info);
		cart_debug("node: %s", uname_info.nodename );
	}
	#endif

	/* set up mpi datatypes, timers, units, etc 
	 * (init_units called in case units set in config file) */
	init_auxiliary();
	init_units();
	init_timers();
	init_logging( restart );
	init_cell_buffer();

#ifdef PARTICLES
	init_particles();
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	init_hydro_tracers();
#endif /* HYDRO_TRACERS */

#ifdef RADIATIVE_TRANSFER
	rtInitRun();
#endif /* RADIATIVE_TRANSFER */

	start_time( TOTAL_TIME );
	start_time( INIT_TIMER );

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
#else
		abox[min_level] = auni[min_level];
#endif
		dt = dtl[min_level];

		for ( i = min_level+1; i <= max_level; i++ ) {
			tl[i] = tl[min_level];
			dtl[i] = 0.5*dtl[i-1];
			abox[i] = abox[min_level];
			auni[i] = auni[min_level];
		}

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
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */

#ifdef GRAVITY
#ifdef HYDRO
		for ( level = min_level; level <= max_level; level++ ) {
			copy_potential( level );
		}
#endif /* HYDRO */
#endif /* GRAVITY */

		check_map();
	} else {
		read_restart(restart_a);
		load_balance(); 
		check_map();

		choose_timestep( &dtl[min_level] );
#ifdef COSMOLOGY
		abox[min_level] = abox_from_tcode(tl[min_level]);
		auni[min_level] = auni_from_tcode(tl[min_level]);
#else
		abox[min_level] = auni[min_level];
#endif
		dt = dtl[min_level];

		for ( i = min_level+1; i <= max_level; i++ ) {
			tl[i] = tl[min_level];
			dtl[i] = 0.5*dtl[i-1];
			abox[i] = abox[min_level];
			auni[i] = auni[min_level];
		}

		for ( i = min_level; i <= max_level; i++ ) {
			cart_debug("num_cells_per_level[%u] = %u", i, num_cells_per_level[i] );
		}

#if defined(RADIATIVE_TRANSFER)
		for ( level = min_level; level <= max_level; level++ ) {
			cart_debug("assigning density on level %u", level );
			assign_density( level );
		}
#endif /* defined(RADIATIVE_TRANSFER) */
	} 

#ifdef COOLING
	init_cooling();
#endif /* COOLING */

	end_time( INIT_TIMER );

	current_steps = 0;
	last_restart_step = 0;

	run_output();

	while ( 1 ) {

#ifdef COSMOLOGY
		cart_debug("a = %e, dt = %e", auni[min_level], dt  );

		if ( auni[min_level] >= auni_end ) {
			break;
		}
#else
		cart_debug("t = %e", tl[min_level] );

		if ( tl[min_level] >= t_end ) {
			break;
		}
#endif /* COSMOLOGY */

		if ( ( timelimit > 0.0 && current_steps > 0 && 
				timelimit-current_time(TOTAL_TIME,min_level) < 1.5*last_time( LEVEL_TIMER, min_level ) ) ||
				( max_steps > 0 && current_steps >= max_steps ) ) {
			cart_debug("reached time or step limit... writing restart");

			/* avoid writing restart file twice on last step */
			if ( last_restart_step != step ) {
				start_time( RESTART_TIMER );
				destroy_cell_buffer();
				write_restart( WRITE_GENERIC, WRITE_GENERIC, WRITE_GENERIC );
				end_time( RESTART_TIMER );
			}

			/* do requeue command */
			if ( local_proc_id == MASTER_NODE && 
					strlen(requeue_command) > 0 ) {
				system( requeue_command );
			}
			break;
		}

		if ( global_timestep( dt ) == -1 ) {
			cart_debug("Error: could not complete timestep, restarting from previous timestep" );

			restart_dt = dtl[min_level];
			choose_timestep( &restart_dt );

			start_time( RESTART_TIMER );
			destroy_cell_buffer();
#ifdef PARTICLES
			init_particles();
#endif /* PARTICLES */
			read_restart(0);
			load_balance(); 
			end_time( RESTART_TIMER );

			dt = 0.5*min( dtl[min_level], restart_dt );
		
#if defined(RADIATIVE_TRANSFER)
			for ( level = min_level; level <= max_level; level++ ) {
			  cart_debug("assigning density on level %u", level );
			  assign_density( level );
			}
#endif /* defined(RADIATIVE_TRANSFER) */
	
		} else {
			step++;
			current_steps++;

			cart_debug("done with timestep %u at tl = %e", step, tl[min_level] );
                        log_diagnostics();

			if ( output_frequency > 0 && step % output_frequency == 0 ) {
				start_time( OUTPUT_TIMER );
				run_output();
				end_time( OUTPUT_TIMER );
			}

			save_check();

			choose_timestep( &dt );
		
			if ( load_balance_frequency > 0 && step % load_balance_frequency == 0 ) {
				load_balance();
			}
		}
	}

	/* destroy buffer */
	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}
	
	cart_debug("total time = %f", end_time(TOTAL_TIME) );

	finalize_logging();
	MPI_Finalize();

	return 0;
}
