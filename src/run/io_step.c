#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include "io.h"
#include "io_cart.h"
#include "io_artio.h"
#include "parallel.h"
#include "rand.h"
#include "times.h"
#include "timing.h"

#include "step.h"

#ifdef LOG_STAR_CREATION
#include "logging.h"
#endif


extern int old_cart_io_flag;
extern int restart_frequency;
extern int grid_output_frequency;

#ifdef PARTICLES
extern int particle_output_frequency;
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
extern int tracer_output_frequency;
#endif /* HYDRO_TRACERS */


void write_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag ) {
	char filename[256];

#ifdef LOG_STAR_CREATION
	char filename_sclog[256];
#endif

	start_time( IO_TIMER );

	if ( old_cart_io_flag ) {
		write_cart_restart( grid_filename_flag, particle_filename_flag, tracer_filename_flag );
	} else {
		write_artio_restart( grid_filename_flag, particle_filename_flag, tracer_filename_flag );
	}

#ifdef STAR_FORMATION
#ifdef LOG_STAR_CREATION
	log_star_creation(-1,-1.0,FILE_CLOSE); //close temp_star files 
	if ( local_proc_id == MASTER_NODE && grid_filename_flag != NO_WRITE ) { 
	  switch(grid_filename_flag) //same switch as hydro which is required anyway for LOG_STAR_CREATION
	    {
	    case WRITE_SAVE:
	      {
#ifdef COSMOLOGY
		sprintf( filename_sclog, "%s/%s_a%06.4f.dsc", output_directory, jobname, auni[min_level] );
#else
		sprintf( filename_sclog, "%s/%s_a%06d.dsc", output_directory, jobname, step );
		cart_error("LOG_STAR_CREATION isn't setup to run without cosmology yet");
#endif /* COSMOLOGY */
		combine_star_creation_log(); 
		finalize_star_creation_log( filename_sclog );
		break;
	      }
	    case WRITE_BACKUP:
	      {
		combine_star_creation_log(); 
		break;
	      }
	    }
	}
	MPI_Barrier( mpi.comm.run ) ; //prevents files from reopening before append occurs
	log_star_creation(-1,-1.0,FILE_OPEN); //close temp_star files 
#endif /* LOG_STAR_CREATION */
#endif /* STAR_FORMATION */

	if ( grid_filename_flag != NO_WRITE && particle_filename_flag != NO_WRITE && tracer_filename_flag != NO_WRITE ) {
		sprintf( filename, "%s/rng_state_"ART_PROC_FORMAT".dat", logfile_directory, local_proc_id );
		cart_rand_save_state( filename );

		last_restart_step = step;
	}

	end_time ( IO_TIMER );
}

void save_check() {
	int grid_save_flag;
	int particle_save_flag;
	int tracer_save_flag;

	grid_save_flag = particle_save_flag = tracer_save_flag = NO_WRITE;

#ifdef COSMOLOGY
	if ( current_output < num_outputs && auni[min_level] >= outputs[current_output] ) {
		grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_SAVE;
		current_output++;
	} else
#endif /* COSMOLOGY */
	if ( restart_frequency != 0 && step % restart_frequency == 0 ) {
		grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_BACKUP;
	} 

	if ( grid_output_frequency != 0 && step % grid_output_frequency == 0 ) {
		grid_save_flag = WRITE_SAVE;		
    }

#ifdef PARTICLES
	if ( particle_output_frequency != 0 && step % particle_output_frequency == 0 ) {
		particle_save_flag = WRITE_SAVE;
	} 
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	if ( tracer_output_frequency != 0 && step % tracer_output_frequency == 0 ) {
		tracer_save_flag = WRITE_SAVE;
	}
#endif /* HYDRO_TRACERS */

	write_restart( grid_save_flag, particle_save_flag, tracer_save_flag );
}
