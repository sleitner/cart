#include "config.h"

#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../base/auxiliary.h"

#include "../core/agn.h"
#include "../core/cell_buffer.h"
#include "../core/cooling.h"
#include "../core/hydro_tracer.h"
#include "../core/io.h"
#include "../core/parallel.h"
#include "../core/particle.h"
#include "../core/rand.h"
#include "../core/rt.h"
#include "../core/starformation.h"
#include "../core/timing.h"

#include "compute_engine.h"


/*
//  Reading control parameters from the config file
*/
void config_init();
void config_read_file(const char *filename);
void config_create_file(const char *filename);
void config_print_to_file(const char *filename, int append);

int drive_run();
int drive_fft();


int drive() {
	int ret;

	configure_runtime_setup();

	if(mpi.task_type & MPI_TASK_TYPE_RUN)
	  {
	    ret = drive_run();
	  }
	else if(mpi.task_type & MPI_TASK_TYPE_FFT)
	  {
	    ret = drive_fft();
	  }
	else
	  {
	    cart_debug("Hooray! I am on vacation this time!");
	    ret = 0;
	  }

	cart_free(options);

	MPI_Finalize();

	return ret;
}


int drive_run () {
	int i;
	int restart = 0;
	char c, *restart_label;
	double aexp;
	const char *str;

        MPI_Comm_size( mpi.comm.run, &num_procs );
        MPI_Comm_rank( mpi.comm.run, &local_proc_id );
        
        if ( num_procs > MAX_PROCS ) {
                cart_error("Number of processors exceeds limit! (%u > %u) The executable must be recompiled.\n", 
                        num_procs, MAX_PROCS );
        }

	/* load configuration file */

	config_init();

	/*
	//  Look for the config file name
	*/
	for(i=0; i<num_options; i++)
	  {
	    if(options[i][0] != '-') break;
	  }
	if(i == num_options)
	  {
	    config_create_file("sample.cfg");
	    cart_error("Usage: art <config_file> [command-line-options]\n   A documented sample of <config_file> is created\n   in the current working directory as sample.cfg");
	  }

	config_read_file(options[i]);
	config_print_to_file("config.log",0);

	/*
	//  Eat out the config file name, we already used it.
	*/
	num_options--;
	for(; i<num_options; i++) options[i] = options[i+1];

	/*
	//  Also support an option-style restart in the form
	//    -r/--restart[=<value>]
	//  where <value> is the value of the scale factor that
	//  labels restart files.
	*/
	str = extract_option1("restart","r","last");
	if(str != NULL)
	  {
	    if(strcmp(str,"last") == 0)
	      {
		restart = 1;
		restart_label = NULL;
	      }
	    else
	      {
		restart = 2;
		restart_label = strstr(str,str);
#ifdef COSMOLOGY
		if(sscanf(str,"%lg%c",&aexp,&c)==1 && aexp>0.0 && aexp<1.1)
		  {
		    /*
		    //  This is a scale factor
		    */
		    restart_label = cart_alloc(char,strlen(str)+1);
		    strcpy(restart_label+1,str);
		    restart_label[0] = 'a';
		  }
#endif
	      }
	  }

	/* start actual work */
	init_rand();
	init_timers();

	start_time( TOTAL_TIME );
	start_time( INIT_TIMER );

	/* set up mpi datatypes, timers, units, etc */
	init_cell_buffer();

#ifdef PARTICLES
	init_particles();
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	init_hydro_tracers();
#endif /* HYDRO_TRACERS */

#ifdef STARFORM
	init_star_formation();
#ifdef AGN
	init_agn();
#endif /* AGN */
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER
	rtInitRun();
#else
#ifdef COOLING
	init_cooling();
#endif /* COOLING */
#endif /* RADIATIVE_TRANSFER */

	run(restart,restart_label);

	return 0;
}


int drive_fft () {
	cart_error("This function is not implemented");
	return -1;
}