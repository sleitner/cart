#include "config.h"

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro_tracer.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_io.h"
#include "sfc.h"
#include "starformation.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


#include "run/hydro_step.h"
#include "run/step.h"

#ifdef LOG_STAR_CREATION
#include "run/logging.h"
#endif


char output_directory[CONTROL_PARAMETER_STRING_LENGTH];
char logfile_directory[CONTROL_PARAMETER_STRING_LENGTH];
char jobname[CONTROL_PARAMETER_STRING_LENGTH];
char requeue_command[CONTROL_PARAMETER_STRING_LENGTH];

int current_output;
int last_restart_step;

int output_frequency = 0;
int restart_frequency = 1;
int particle_output_frequency = 0;
int tracer_output_frequency = 0;
int grid_output_frequency = 0;

int num_outputs = 0;
int outputs_size = 0;
float *outputs = NULL;

int old_art_io_flag = 0;

/*
//  For backward compatibility
*/
void units_set_art(double OmegaM, double h, double Lbox);


void control_parameter_set_outputs(const char *value, void *ptr, int ind)
{
  int i, j;
  char *str, *tok;
  float a1, a2, da, *tmp;

  str = cart_alloc(char,strlen(value)+1);
  strcpy(str,value); /* strtok destroys the input string */

  tok = strtok(str," ");
  while(tok != NULL)
    {
      if(sscanf(tok,"(%g,%g,%g)",&a1,&a2,&da) == 3)
	{
	  cart_assert(da > 0.0);
	  while(a1 < a2)
	    {
	      if(num_outputs == outputs_size)
		{
		  outputs_size += 100;
		  tmp = cart_alloc(float,outputs_size);
		  memcpy(tmp,outputs,num_outputs*sizeof(float));
		  cart_free(outputs);
		  outputs = tmp;
		}
	      outputs[num_outputs++] = a1;
	      a1 += da;
	    }
	}
      else if(sscanf(tok,"%g",&a1) == 1)
	{
	  if(num_outputs == outputs_size)
	    {
	      outputs_size += 100;
	      tmp = cart_alloc(float,outputs_size);
	      memcpy(tmp,outputs,num_outputs*sizeof(float));
	      cart_free(outputs);
	      outputs = tmp;
	    }
	  outputs[num_outputs++] = a1;
	}
      else
	{
	  cart_error("Invalid snapshot token %s in line %s",tok,value);
	}
      tok = strtok(NULL," ");
    }

  cart_free(str);

  qsort(outputs,num_outputs,sizeof(float),compare_floats);
  /*
  //  Remove duplicates
  */
  for(i=1; i<num_outputs; i++)
    {
      if(outputs[i] < outputs[i-1]+1.0e-6)
	{
	  num_outputs--;
	  for(j=i; j<num_outputs; j++) 
	    {
	      outputs[j] = outputs[j+1];
	    }
	}
    }
}


void control_parameter_print_name(FILE *f, const char *name);

void control_parameter_list_outputs(FILE *stream, const void *ptr)
{
  const int num_per_line = 10;
  int i;
  int newline = 0;

  for(i=0; i<num_outputs; i++)
    {
      if(newline)
	{
	  fprintf(stream,"\n");
	  control_parameter_print_name(stream,"snapshot-epochs");
	  newline = 0;
	}

      fprintf(stream,"%g ",outputs[i]);

      if((i+1)%num_per_line == 0)  newline = 1;
    }

  if(num_outputs == 0)
    {
      /*
      //  Refinement indicators are not set
      */
      //fprintf(stream,"(NOT SET)");
    }
}


void config_init_io()
{
  ControlParameterOps control_parameter_outputs = { control_parameter_set_outputs, control_parameter_list_outputs };

  outputs_size = 100;
  outputs = cart_alloc(float,outputs_size);

  strcpy(jobname,"ART");
  strcpy(output_directory,".");
  strcpy(logfile_directory,".");
  strcpy(requeue_command,"");

  control_parameter_add2(control_parameter_string,jobname,"job-name","jobname","a name for this job. This name can be stored in output files for further reference.");

  control_parameter_add2(control_parameter_string,output_directory,"directory:outputs","output_directory","directory for output files.");

  control_parameter_add2(control_parameter_string,logfile_directory,"directory:logs","logfile_directory","directory for logs and other auxiliary files.");

  control_parameter_add2(control_parameter_string,requeue_command,"requeue-command","requeue_command","a command for re-queueing the batch job after the completion of the current one.");

#ifdef COSMOLOGY
  control_parameter_add2(control_parameter_outputs,outputs,"snapshot-epochs","outputs","values for the cosmic scale factor at which to produce the full snalshots from the simulation. Values can be listed either one by one, separated by white space, or in a form '(a1,a2,da)' which expand into a loop from a=a1 to a=a2, with the step da. For example, '(0.1,0.5,0.1)' is equivalent to '0.1 0.2 0.3 0.4 0.5'. An arbitrary number of entries is allowed, so the loop entry does not have to be just one, several of them can be combined together; for example, '(0.1,0.3,0.1) (0.4,1.0,0.2)' will expand to '0.1 0.2 0.3 0.4 0.6 0.8 1.0'. Entries are sorted and repeated entries are automatically deleted, so the user does not need to care about the order of entries or the existence of duplicates.");
#endif /* COSMOLOGY */

  control_parameter_add2(control_parameter_int,&output_frequency,"frequency:user-output","output_frequency","frequency (in global time steps) of calling a user-defined run_output() function. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&restart_frequency,"frequency:restart","restart_frequency","frequency (in global time steps) of producing restart files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&particle_output_frequency,"frequency:particle-output","particle_output_frequency","frequency (in global time steps) of producing particle output files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&tracer_output_frequency,"frequency:tracer-output","tracer_output_frequency","frequency (in global time steps) of producing tracer output files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&grid_output_frequency,"frequency:grid-output","grid_output_frequency","frequency (in global time steps) of calling producing grid output files. Every time a grid file is writtent to disk, particle and tracer files are also written, irrespectively of the values of <frequency:particle-output> and <frequency:tracer-output> parameters. Zero frequency disables regular output of grid files; however, outputs are still produced in cosmological simulations at cosmic scale factors set by <snapshot-epochs> parameter.");

  control_parameter_add2(control_parameter_bool,&old_art_io_flag,"io:use-old-art-format","old-art-format","Disable using the new (default) I/O library in favor of the original ART particle and grid formats (either T or F).");

  config_init_io_art();
  config_init_io_cartio();
}


void config_verify_io()
{
  int i;
  DIR *d;
  
  if((d = opendir(output_directory)) == NULL)
    {
      cart_error("Directory %s does not exist.",output_directory);
    }
  else closedir(d);

  if((d = opendir(logfile_directory)) == NULL)
    {
      cart_error("Directory %s does not exist.",logfile_directory);
    }
  else closedir(d);

#ifdef COSMOLOGY
  for(i=1; i<num_outputs; i++) if(!(outputs[i-1] < outputs[i]))
    {
      cart_error("Outputs are not strictly increasing (%d,%f,%f)",i,outputs[i-1],outputs[i]);
    }
#endif /* COSMOLOGY */

  cart_assert(output_frequency >= 0);

  cart_assert(restart_frequency >= 0);

  cart_assert(particle_output_frequency >= 0);

  cart_assert(tracer_output_frequency >= 0);

  cart_assert(grid_output_frequency >= 0);

  cart_assert(old_art_io_flag == 0 || old_art_io_flag == 1);

  config_verify_io_art();
  config_verify_io_cartio();
}

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
		save_auxiliary();
		last_restart_step = step;
	}

	end_time ( IO_TIMER );
}

void read_restart( const char *label ) {
	start_time( IO_TIMER );

	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}

	if ( old_art_io_flag ) {
		read_art_restart(label);	
	} else {
		read_cartio_restart(label);
	}

	units_reset();
	units_update(min_level);

#ifdef COSMOLOGY
	/* ensure current output points to correct place in output array */
	current_output = 0;
	while ( current_output < num_outputs &&
			auni[min_level] >= outputs[current_output] ) {
		current_output++;
	}
#endif /* COSMOLOGY */

	end_time( IO_TIMER );
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
