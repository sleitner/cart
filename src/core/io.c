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
#include "io_cart.h"
#include "io_artio.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "rand.h"
#include "refinement.h"
#include "rt_io.h"
#include "sfc.h"
#include "starformation.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


extern int step;


char output_directory_d[CONTROL_PARAMETER_STRING_LENGTH] = ".";
const char* output_directory = output_directory_d;

char logfile_directory_d[CONTROL_PARAMETER_STRING_LENGTH] = ".";

char jobname_d[CONTROL_PARAMETER_STRING_LENGTH] = "ART";
const char* jobname = jobname_d;

char requeue_command_d[CONTROL_PARAMETER_STRING_LENGTH] = "";
const char* requeue_command = requeue_command_d;

int current_output;
int last_restart_step;

int current_restart_backup = 0;
int num_restart_backups = 2;

int output_frequency = 0;
int restart_frequency = 1;
#ifdef PARTICLES
int particle_output_frequency = 0;
#endif /* PARTICLES */
#ifdef HYDRO_TRACERS
int tracer_output_frequency = 0;
#endif /* HYDRO_TRACERS */
int grid_output_frequency = 0;

int num_outputs = 0;
int outputs_size = 0;
float *outputs = NULL;

int old_cart_io_flag = 0;

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
	  while(a1 < a2+0.5*da)
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

  control_parameter_add2(control_parameter_string,jobname_d,"job-name","jobname","a name for this job. This name can be stored in output files for further reference.");

  control_parameter_add2(control_parameter_string,output_directory_d,"directory:outputs","output_directory","directory for output files.");

  control_parameter_add2(control_parameter_string,logfile_directory_d,"directory:logs","logfile_directory","directory for logs and other auxiliary files.");

  control_parameter_add2(control_parameter_string,requeue_command_d,"requeue-command","requeue_command","a command for re-queueing the batch job after the completion of the current one.");

#ifdef COSMOLOGY
  control_parameter_add2(control_parameter_outputs,outputs,"snapshot-epochs","outputs","values for the cosmic scale factor at which to produce the full snalshots from the simulation. Values can be listed either one by one, separated by white space, or in a form '(a1,a2,da)' which expand into a loop from a=a1 to a=a2, with the step da. For example, '(0.1,0.5,0.1)' is equivalent to '0.1 0.2 0.3 0.4 0.5'. An arbitrary number of entries is allowed, so the loop entry does not have to be just one, several of them can be combined together; for example, '(0.1,0.3,0.1) (0.4,1.0,0.2)' will expand to '0.1 0.2 0.3 0.4 0.6 0.8 1.0'. Entries are sorted and repeated entries are automatically deleted, so the user does not need to care about the order of entries or the existence of duplicates.");
#endif /* COSMOLOGY */

  control_parameter_add2(control_parameter_int,&output_frequency,"frequency:user-output","output_frequency","frequency (in global time steps) of calling a user-defined run_output() function. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&restart_frequency,"frequency:restart","restart_frequency","frequency (in global time steps) of producing restart files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&num_restart_backups,"io:restart-backups","num_restart_backups","number of individual restart backups to retain before overwriting.");

#ifdef PARTICLES
  control_parameter_add2(control_parameter_int,&particle_output_frequency,"frequency:particle-output","particle_output_frequency","frequency (in global time steps) of producing particle output files. Zero frequency disables this option.");
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
  control_parameter_add2(control_parameter_int,&tracer_output_frequency,"frequency:tracer-output","tracer_output_frequency","frequency (in global time steps) of producing tracer output files. Zero frequency disables this option.");
#endif /* HYDRO_TRACERS */

  control_parameter_add2(control_parameter_int,&grid_output_frequency,"frequency:grid-output","grid_output_frequency","frequency (in global time steps) of calling producing grid output files. Every time a grid file is writtent to disk, particle and tracer files are also written, irrespectively of the values of <frequency:particle-output> and <frequency:tracer-output> parameters. Zero frequency disables regular output of grid files; however, outputs are still produced in cosmological simulations at cosmic scale factors set by <snapshot-epochs> parameter.");

  control_parameter_add2(control_parameter_bool,&old_cart_io_flag,"io:use-old-cart-format","old-cart-format","Disable using the new (default) I/O library in favor of the original ART particle and grid formats (either T or F).");

  config_init_io_cart();
  config_init_io_artio();
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

  VERIFY(snapshot-epochs, 1 );
#endif /* COSMOLOGY */

  VERIFY(frequency:user-output, output_frequency >= 0 );

  VERIFY(frequency:restart, restart_frequency >= 0 );

  VERIFY(io:restart-backups, num_restart_backups > 0 );

#ifdef PARTICLES
  VERIFY(frequency:particle-output, particle_output_frequency >= 0 );
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
  VERIFY(frequency:tracer-output, tracer_output_frequency >= 0 );
#endif /* HYDRO_TRACERS */

  VERIFY(frequency:grid-output, grid_output_frequency >= 0 );

  config_verify_io_cart();
  config_verify_io_artio();
}

void read_restart( const char *label ) {
	char filename[256];

	start_time( IO_TIMER );

	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}

	if ( old_cart_io_flag ) {
		read_cart_restart(label);	
		if ( old_cart_io_flag < 0 ) old_cart_io_flag = 0;
	} else {
		read_artio_restart(label);
	}

	/* load random number generator state */
	sprintf( filename, "%s/rng_state_"ART_PROC_FORMAT".dat", logfile_directory, local_proc_id );
	cart_rand_load_state( filename, 0 );

	units_init();
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


void set_jobname(const char *str)
{
  if(str != NULL)
    {
      strncpy(jobname_d,str,CONTROL_PARAMETER_STRING_LENGTH);
      jobname_d[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
    }
}
