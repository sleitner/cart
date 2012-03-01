#include "config.h"


#include "io.h"
#include "io_art.h"
#include "io_cartio.h"
#include "rand.h"
#include "times.h"
#include "timing.h"

#include "step.h"

#ifdef LOG_STAR_CREATION
#include "logging.h"
#endif


extern int old_art_io_flag;
extern int restart_frequency;
extern int grid_output_frequency;
extern int particle_output_frequency;
extern int tracer_output_frequency;


void write_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag ) {
#ifdef LOG_STAR_CREATION
	char filename_sclog[256];
#endif

	start_time( IO_TIMER );

	if ( old_art_io_flag ) {
		write_art_restart( grid_filename_flag, particle_filename_flag, tracer_filename_flag );
	} else {
		write_cartio_restart( grid_filename_flag, particle_filename_flag, tracer_filename_flag );
	}

#ifdef STARFORM
#ifdef LOG_STAR_CREATION
	log_star_creation(-1,-1.0,FILE_CLOSE); //close temp_star files 
	if ( local_proc_id == MASTER_NODE && grid_filename_flag != NO_WRITE ) { 
	  switch(grid_filename_flag) //same switch as hydro which is required anyway for LOG_STAR_CREATION
	    {
	    case WRITE_SAVE:
	      {
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
#ifdef COSMOLOGY
		sprintf( filename_sclog, "%s/%s_a%06.4f.dsc", output_directory, jobname, auni[min_level] );
#else
		sprintf( filename_sclog, "%s/%s_a%06d.dsc", output_directory, jobname, step );
		cart_error("LOG_STAR_CREATION isn't setup to run without cosmology yet");
#endif /* COSMOLOGY */
#else
#ifdef COSMOLOGY
		sprintf( filename_sclog, "%s/star_creation_a%06.4f.dat", output_directory, auni[min_level] );
#else
		sprintf( filename_sclog, "%s/star_creation_a%06d.dat", output_directory, step );
		cart_error("LOG_STAR_CREATION define isn't setup to run without cosmology yet");
#endif /* COSMOLOGY */
#endif
		combine_star_creation_log(); 
		finalize_star_creation_log( filename_sclog );
		break;
	      }
	    case WRITE_BACKUP:
	      {
		//sprintf( filename_sclog, "%s/%s_2.dsc", output_directory, jobname );
		combine_star_creation_log(); 
		break;
	      }
	    case WRITE_GENERIC:
	      {
		//sprintf( filename_sclog, "%s/%s.dsc", output_directory, jobname );
		combine_star_creation_log(); 
		break;
	      }
	    }
	}
	MPI_Barrier( mpi.comm.run ) ; //prevents files from reopening before append occurs
	log_star_creation(-1,-1.0,FILE_OPEN); //close temp_star files 
#endif /* LOG_STAR_CREATION */
#endif /* STARFORM */

	if ( grid_filename_flag != NO_WRITE ) {
		save_rand();
		last_restart_step = step;
	}

	end_time ( IO_TIMER );
}

void save_check() {
	int grid_save_flag;
	int particle_save_flag;
	int tracer_save_flag;

	grid_save_flag = particle_save_flag = tracer_save_flag = NO_WRITE;

	if ( restart_frequency != 0 && step % restart_frequency == 0 ) {
		if ( step % (2*restart_frequency) == 0 ) {
			grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_BACKUP;
		} else {
			grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_GENERIC;
		}
#ifdef COSMOLOGY
	} else if ( current_output < num_outputs && auni[min_level] >= outputs[current_output] ) {
		grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_SAVE;
		current_output++;
#endif /* COSMOLOGY */
	} 

	
	if ( grid_output_frequency != 0 && step % grid_output_frequency == 0 ) {
		grid_save_flag = WRITE_SAVE;		
    }
 
	if ( particle_output_frequency != 0 && step % particle_output_frequency == 0 ) {
		particle_save_flag = WRITE_SAVE;
	} 

	if ( tracer_output_frequency != 0 && step % tracer_output_frequency == 0 ) {
		tracer_save_flag = WRITE_SAVE;
	}

	write_restart( grid_save_flag, particle_save_flag, tracer_save_flag );
}
