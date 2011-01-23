#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cooling.h"
#include "cosmology.h"
#include "density.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "io.h"
#include "load_balance.h"
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "plugin.h"
#include "rt_solver.h"
#include "starformation.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"
#include "top_level_fft.h"


/*
//  Reading control parameters from the config file
*/
void config_init();
void config_read_file(const char *filename);
void config_create_file(const char *filename);
void config_print_to_file(const char *filename, int append);
void config_append_units_to_file(const char *filename);


void init_run();
void run_output();
void config_plugins();


int main ( int argc, char *argv[]) {
	int current_steps;
	int restart;
	int level, mode;
	char c, *restart_label;
	double aexp;
	const char *tmp;
#ifdef _OPENMP
	int nomp;
#endif

	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &num_procs );
	MPI_Comm_rank( MPI_COMM_WORLD, &local_proc_id );
	
	if ( num_procs > MAX_PROCS ) {
		cart_error("Number of processors exceeds limit! (%u > %u) The executable must be recompiled.\n", 
			num_procs, MAX_PROCS );
	}

	/* load configuration file */
	init_auxiliary();
	init_timers();

	start_time( TOTAL_TIME );
	start_time( INIT_TIMER );

	config_init();
	if ( argc < 2 )
	  {
	    config_create_file("sample.cfg");
	    cart_error("Usage: art <config_file> [-r/--restart[=<aexp>]] [-pfm/--particle-file-mode=<mode>] [-gfm/--grid-file-mode=<mode>] ...\n   A documented sample of <config_file> is created\n   in the current working directory as sample.cfg");
	  }
	else
	  {
	    config_read_file( argv[1] );
	  }

	config_print_to_file("config.log",0);

	if ( argc == 2 ) {
		restart = 0;
		/* skip config file name */
		num_options = argc - 2;
		options = argv + 2;
	} else {
		/*
		//  Also support an option-style restart in the form
		//    -r/--restart[=<value>]
		//  where <value> is the value of the scale factor that
		//  labels restart files.
		*/
		tmp = check_option1(argv[2],"-restart","last");
		if(tmp == NULL) tmp = check_option1(argv[2],"r","last");
		if(tmp == NULL)
		  {
		    restart = 0;
		    /* skip config file name */
		    num_options = argc - 2;
		    options = argv + 2;
		  }
		else
		  {
		    if(strcmp(tmp,"last") == 0)
		      {
			restart = 1;
			restart_label = NULL;
		      }
		    else
		      {
			restart = 2;
			restart_label = strstr(argv[2],tmp);
#ifdef COSMOLOGY
			if(sscanf(tmp,"%lg%c",&aexp,&c)==1 && aexp>0.0 && aexp<1.1)
			  {
			    /*
			    //  This is a scale factor
			    */
			    restart_label = cart_alloc(char,strlen(tmp)+1);
			    strcpy(restart_label+1,tmp);
			    restart_label[0] = 'a';
			  }
#endif
		      }
		    /* skip config file name and restart flag */
		    num_options = argc - 3;
		    options = argv + 3;
		  }

#ifdef _OPENMP
		/*
		//  Support -omp/--num-omp-threads=<num-threads> option
		*/
		if(num_options > 0)
		  {
		    tmp = check_option1(options[0],"-num-omp-threads",NULL);
		    if(tmp == NULL) tmp = check_option1(options[0],"omp",NULL);
		    if(tmp != NULL)
		      {
			if(sscanf(tmp,"%d",&nomp)!=1 || nomp<1 || nomp>256)
			  {
			    cart_error("-omp=<num> option requires a positive integer <num> as an argument");
			  }
			omp_set_num_threads(nomp);
			num_options--;
			options++;
		      }
		  }
#endif
#ifdef PARTICLES
		/*
		//  Support -pfm/--particle-file-mode=<mode> option
		*/
		if(num_options > 0)
		  {
		    tmp = check_option1(options[0],"-particle-file-mode",NULL);
		    if(tmp == NULL) tmp = check_option1(options[0],"pfm",NULL);
		    if(tmp != NULL)
		      {
			if(sscanf(tmp,"%d",&mode)==0 || mode<0 || mode>2)
			  {
			    cart_error("-pfm=<mode> requires an integer <mode> between 0 and 2 as an argument.");
			  }
			set_particle_file_mode(mode);
			options++;
			num_options--;
		      }
		  }
#endif /* PARTICLES */
		/*
		//  Support -gfm/--grid-file-mode=<mode> option
		*/
		if(num_options > 0)
		  {
		    tmp = check_option1(options[0],"-grid-file-mode",NULL);
		    if(tmp == NULL) tmp = check_option1(options[0],"gfm",NULL);
		    if(tmp != NULL)
		      {
			if(sscanf(tmp,"%d",&mode)==0 || mode<0 || mode>3)
			  {
			    cart_error("-gfm=<mode> requires an integer <mode> between 0 and 3 as an argument.");
			  }
			set_grid_file_mode(mode);
			options++;
			num_options--;
		      }
		  }
	}

#ifdef _OPENMP
	cart_debug("num openmp threads = %u", omp_get_max_threads() );
#endif

	/* set up mpi datatypes, timers, units, etc */
	init_logging( restart );
	init_cell_buffer();

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER) 
	init_fft();
#endif

#ifdef PARTICLES
	init_particles();
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	init_hydro_tracers();
#endif /* HYDRO_TRACERS */

#ifdef STARFORM
	init_star_formation();
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER
	rtInitRun();
#else
#ifdef COOLING
	init_cooling();
#endif /* COOLING */
#endif /* RADIATIVE_TRANSFER */

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

	end_time( INIT_TIMER );

	current_steps = 0;
	last_restart_step = 0;

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

			read_restart(0);
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
	MPI_Finalize();

	return 0;
}
