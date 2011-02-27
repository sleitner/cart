#include "config.h"

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "index_hash.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "logging.h" 
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_io.h"
#include "sfc.h"
#include "skiplist.h"
#include "starformation.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


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

int num_output_files = 1;

int num_outputs = 0;
int outputs_size = 0;
float *outputs = NULL;


#ifndef PARTICLE_HEADER_MAGIC
#define PARTICLE_HEADER_MAGIC           (0.1234f)
#endif


int grid_file_mode = 0;
int particle_file_mode = 0;


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

  control_parameter_add2(control_parameter_int,&num_output_files,"num-output-files","num_output_files","Number of parallel output files.");

#ifdef COSMOLOGY
  control_parameter_add2(control_parameter_outputs,outputs,"snapshot-epochs","outputs","values for the cosmic scale factor at which to produce the full snalshots from the simulation. Values can be listed either one by one, separated by white space, or in a form '(a1,a2,da)' which expand into a loop from a=a1 to a=a2, with the step da. For example, '(0.1,0.5,0.1)' is equivalent to '0.1 0.2 0.3 0.4 0.5'. An arbitrary number of entries is allowed, so the loop entry does not have to be just one, several of them can be combined together; for example, '(0.1,0.3,0.1) (0.4,1.0,0.2)' will expand to '0.1 0.2 0.3 0.4 0.6 0.8 1.0'. Entries are sorted and repeated entries are automatically deleted, so the user does not need to care about the order of entries or the existence of duplicates.");
#endif /* COSMOLOGY */

  control_parameter_add2(control_parameter_int,&output_frequency,"frequency:user-output","output_frequency","frequency (in global time steps) of calling a user-defined run_output() function. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&restart_frequency,"frequency:restart","restart_frequency","frequency (in global time steps) of producing restart files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&particle_output_frequency,"frequency:particle-output","particle_output_frequency","frequency (in global time steps) of producing particle output files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&tracer_output_frequency,"frequency:tracer-output","tracer_output_frequency","frequency (in global time steps) of producing tracer output files. Zero frequency disables this option.");

  control_parameter_add2(control_parameter_int,&grid_output_frequency,"frequency:grid-output","grid_output_frequency","frequency (in global time steps) of calling producing grid output files. Every time a grid file is writtent to disk, particle and tracer files are also written, irrespectively of the values of <frequency:particle-output> and <frequency:tracer-output> parameters. Zero frequency disables regular output of grid files; however, outputs are still produced in cosmological simulations at cosmic scale factors set by <snapshot-epochs> parameter.");
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

  cart_assert(num_output_files > 0);

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
}


void reorder( char *buffer, int size ) {
        int i;
        char tmp;
                                                                                
        for ( i = 0; i < (size/2); i++ ) {
                tmp = buffer[i];
                buffer[i] = buffer[size - i - 1];
                buffer[size - i - 1] = tmp;
        }
}

void write_restart( int gas_filename_flag, int particle_filename_flag, int tracer_filename_flag, char *label ) {
	FILE *restart;
	char filename[256];
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
	char filename_gas[256];
#ifdef LOG_STAR_CREATION
	char filename_sclog[256];
#endif
#ifdef HYDRO
#ifdef HYDRO_TRACERS
	char filename_tracers[256];
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */
	char label_buffer[256];

	start_time( IO_TIMER );

	if(label == NULL)
	  {
	    label = label_buffer;
#ifdef COSMOLOGY
	    sprintf(label,"a%06.4f",auni[min_level]);
#else
	    sprintf(label,"%06d",step);
#endif /* COSMOLOGY */
	  }

	cart_debug("Writing grid restart...");
	switch(gas_filename_flag)
	  {
	  case WRITE_SAVE:
	    {
	      sprintf( filename_gas, "%s/%s_%s.d", output_directory, jobname, label );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
	      sprintf( filename_gas, "%s/%s_2.d", output_directory, jobname );
	      break;
	    }
	  default:
	    {
	      sprintf( filename_gas, "%s/%s.d", output_directory, jobname );
	    }
	}

	start_time( GAS_WRITE_IO_TIMER );
	write_grid_binary( filename_gas );
	end_time( GAS_WRITE_IO_TIMER );

#ifdef HYDRO
#ifdef HYDRO_TRACERS
	cart_debug("Writing hydro tracer restart...");
	switch(tracer_filename_flag)
	  {
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
	  case WRITE_SAVE:
	    {
	      sprintf( filename_tracers, "%s/%s_%s.dtr", output_directory, jobname, label );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
	      sprintf( filename_tracers, "%s/%s_2.dtr", output_directory, jobname );
	      break;
	    }
	  default:
	    {
	      sprintf( filename_tracers, "%s/%s.dtr", output_directory, jobname );
	    }
#else  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
	  case WRITE_SAVE:
	    {
	      sprintf( filename_tracers, "%s/tracers_%s.dat", output_directory, label );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
	      sprintf( filename_tracers, "%s/tracers_2.dat", output_directory );
	      break;
	    }
	  default:
	    {
	      sprintf( filename_tracers, "%s/tracers.dat", output_directory );
	    }
#endif /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
	}

	start_time( PARTICLE_WRITE_IO_TIMER );
	write_hydro_tracers( filename_tracers );
	end_time( PARTICLE_WRITE_IO_TIMER );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#ifdef PARTICLES
	cart_debug("Writing particle restart...");
	switch(particle_filename_flag)
	  {
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
	  case WRITE_SAVE:
	    {
	      sprintf( filename1, "%s/%s_%s.dph", output_directory, jobname, label );
	      sprintf( filename2, "%s/%s_%s.dxv", output_directory, jobname, label );
	      sprintf( filename3, "%s/%s_%s.dpt", output_directory, jobname, label );
	      sprintf( filename4, "%s/%s_%s.dst", output_directory, jobname, label );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
		sprintf( filename1, "%s/%s_2.dph", output_directory, jobname);
		sprintf( filename2, "%s/%s_2.dxv", output_directory, jobname);
		sprintf( filename3, "%s/%s_2.dpt", output_directory, jobname);
		sprintf( filename4, "%s/%s_2.dst", output_directory, jobname);
		break;
	    }
	  default:
	    {
		sprintf( filename1, "%s/%s.dph", output_directory, jobname);
		sprintf( filename2, "%s/%s.dxv", output_directory, jobname);
		sprintf( filename3, "%s/%s.dpt", output_directory, jobname);
		sprintf( filename4, "%s/%s.dst", output_directory, jobname);
	    }
#else  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
	  case WRITE_SAVE:
	    {
	      sprintf( filename1, "%s/PMcrd%s.DAT", output_directory, label );
	      sprintf( filename2, "%s/PMcrs0%s.DAT", output_directory, label );
	      sprintf( filename3, "%s/pt%s.dat", output_directory, label );
	      sprintf( filename4, "%s/stars_%s.dat", output_directory, label );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
		sprintf( filename1, "%s/PMcrd_2.DAT", output_directory );
		sprintf( filename2, "%s/PMcrs_2.DAT", output_directory );
		sprintf( filename3, "%s/pt_2.dat", output_directory );
		sprintf( filename4, "%s/stars_2.dat", output_directory );
		break;
	    }
	  default:
	    {
		sprintf( filename1, "%s/PMcrd.DAT", output_directory );
		sprintf( filename2, "%s/PMcrs.DAT", output_directory );
		sprintf( filename3, "%s/pt.dat", output_directory );
		sprintf( filename4, "%s/stars.dat", output_directory );
	    }
#endif  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
	}
		
	start_time( PARTICLE_WRITE_IO_TIMER );
#ifdef STARFORM
	write_particles( filename1, filename2, filename3, filename4 );
#else
	write_particles( filename1, filename2, filename3, NULL );
#endif /* STARFORM */
	end_time( PARTICLE_WRITE_IO_TIMER );

#endif /* PARTICLES */


#ifdef STARFORM
#ifdef LOG_STAR_CREATION
	log_star_creation(-1,-1.0,FILE_CLOSE); //close temp_star files 
	if ( local_proc_id == MASTER_NODE ) { 
	  switch(gas_filename_flag) //same switch as hydro which is required anyway for LOG_STAR_CREATION
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
	    default:
	      {
		//sprintf( filename_sclog, "%s/%s.dsc", output_directory, jobname );
		combine_star_creation_log(); 
	      }
	    }
	}
	MPI_Barrier( MPI_COMM_WORLD ) ; //prevents files from reopening before append occurs
	log_star_creation(-1,-1.0,FILE_OPEN); //close temp_star files 
#endif /*LOG_STAR_CREATION */
#endif /*STARFORM */


	if ( local_proc_id == MASTER_NODE ) {
		/* write out restart file */
		sprintf( filename, "%s/restart.dat", output_directory );
		restart = fopen( filename, "w" );
		cart_assert(restart != NULL);

		fprintf( restart, "%s\n", filename_gas );

#ifdef HYDRO
#ifdef HYDRO_TRACERS
		fprintf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#ifdef PARTICLES
		fprintf( restart, "%s\n", filename1 );
		fprintf( restart, "%s\n", filename2 );
		fprintf( restart, "%s\n", filename3 );
#ifdef STARFORM
		fprintf( restart, "%s\n", filename4 );
#endif /* STARFORM */
#endif /* PARTICLES */
	
		fclose(restart);
	}

	save_auxiliary();

	last_restart_step = step;

	end_time ( IO_TIMER );
}

void read_restart( const char *label ) {
	FILE *restart;
	char filename[256];
	char filename_gas[256];
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
	char filename_tracers[256];

#ifdef SAVE_LOAD_BALANCE_PARTITION
	FILE *partition;
	int partition_num_grid, partition_num_procs, proc;
#endif

	start_time( IO_TIMER );

	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}

	if ( label == NULL ) {
		/* read filenames from restart file */
		sprintf( filename, "%s/restart.dat", output_directory );
		restart = fopen( filename, "r" );

		if ( restart == NULL ) {
			cart_debug("Unable to locate restart.dat, trying default filenames!");

			/* try generic names */
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
			sprintf( filename_gas, "%s/%s.d", output_directory, jobname );
			sprintf( filename1,  "%s/%s.dph", output_directory, jobname );
			sprintf( filename2, "%s/%s.dxv", output_directory, jobname );
			sprintf( filename3, "%s/%s.dpt", output_directory, jobname );
			sprintf( filename4, "%s/%s.dst", output_directory, jobname );
			sprintf( filename_tracers, "%s/%s.dtr", output_directory, jobname );
#else
			sprintf( filename_gas, "%s/%s.d", output_directory, jobname );
			sprintf( filename1,  "%s/PMcrd.DAT", output_directory );
			sprintf( filename2, "%s/PMcrs.DAT", output_directory );
			sprintf( filename3, "%s/pt.dat", output_directory );
			sprintf( filename4, "%s/stars.dat", output_directory );
			sprintf( filename_tracers, "%s/tracers.dat", output_directory );
#endif
		} else {
			fscanf( restart, "%s\n", filename_gas );
#ifdef HYDRO
#ifdef HYDRO_TRACERS
			fscanf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */
			fscanf( restart, "%s\n", filename1 );
			fscanf( restart, "%s\n", filename2 );
			fscanf( restart, "%s\n", filename3 );
#ifdef STARFORM
			fscanf( restart, "%s\n", filename4 );
#endif /* STARFORM */
			fclose(restart);
		}
	} else {
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
		sprintf( filename_gas, "%s/%s_%s.d", output_directory, jobname, label );
		sprintf( filename1,  "%s/%s_%s.dph", output_directory, jobname, label );
		sprintf( filename2, "%s/%s_%s.dxv", output_directory, jobname, label );
		sprintf( filename3, "%s/%s_%s.dpt", output_directory, jobname, label );
		sprintf( filename4, "%s/%s_%s.dst", output_directory, jobname, label );
		sprintf( filename_tracers, "%s/%s_%s.dtr", output_directory, jobname, label );
#else
		sprintf( filename_gas, "%s/%s_%s.d", output_directory, jobname, label );
		sprintf( filename1,  "%s/PMcrd%s.DAT", output_directory, label );
		sprintf( filename2, "%s/PMcrs0%s.DAT", output_directory, label );
		sprintf( filename3, "%s/pt%s.dat", output_directory, label );
		sprintf( filename4, "%s/stars_%s.dat", output_directory, label );
		sprintf( filename_tracers, "%s/tracers_%s.dat", output_directory, label );
#endif
	}

#ifdef SAVE_LOAD_BALANCE_PARTITION
    /* try to find load-balancing file */
    sprintf( filename, "%s/partition.dat", output_directory );                                                                    
    partition = fopen( filename, "r" );

    if ( partition == NULL ) {
        /* do load balancing */
        restart_load_balance( filename_gas, filename1, filename2 );
    } else {
		fscanf( partition, "%u %u\n", &partition_num_procs, &partition_num_grid );

		if ( partition_num_procs == num_procs && partition_num_grid == num_grid ) {
			for ( proc = 0; proc < num_procs+1; proc++ ) {
				fscanf( partition, "%u\n", &proc_sfc_index[proc] );
			}

			init_tree();
		} else {
			cart_error("Bad values for num_procs or num_grid in partition.dat!");
		}

		fclose( partition );
	}
#else
	/* do load balancing */
	restart_load_balance( filename_gas, filename1, filename2 );
#endif

	cart_debug("Reading grid restart...");
	start_time( GAS_READ_IO_TIMER );
	read_grid_binary( filename_gas );
	end_time( GAS_READ_IO_TIMER );

#ifdef HYDRO
#ifdef HYDRO_TRACERS
	start_time( PARTICLE_READ_IO_TIMER );
	read_hydro_tracers( filename_tracers );
	end_time( PARTICLE_READ_IO_TIMER );
#endif /* HYDRO_TRACERS */

#else
	build_root_cell_buffer();
#endif /* HYDRO */

#ifdef PARTICLES
	cart_debug("Reading particle restart...");
	start_time( PARTICLE_READ_IO_TIMER );
#ifdef STARFORM
	read_particles( filename1, filename2, filename3, filename4, 0, NULL );
#else
	read_particles( filename1, filename2, filename3, NULL, 0, NULL );
#endif /* STARFORM */
	end_time( PARTICLE_READ_IO_TIMER );
#endif /* PARTICLES */

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

	/*
	//  If we started from a labeled file and succeeded, write
	//  a restart file - just in case we encounter a CFL violation
	//  in the first step.
	*/
	if ( label != NULL && local_proc_id == MASTER_NODE ) {
		/* write out restart file */
		sprintf( filename, "%s/restart.dat", output_directory );
		restart = fopen( filename, "w" );
		cart_assert(restart != NULL);

		fprintf( restart, "%s\n", filename_gas );

#ifdef HYDRO
#ifdef HYDRO_TRACERS
		fprintf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#ifdef PARTICLES
		fprintf( restart, "%s\n", filename1 );
		fprintf( restart, "%s\n", filename2 );
		fprintf( restart, "%s\n", filename3 );
#ifdef STARFORM
		fprintf( restart, "%s\n", filename4 );
#endif /* STARFORM */
#endif /* PARTICLES */
	
		fclose(restart);
	}

	end_time( IO_TIMER );
}

void save_check() {
	char filename1[256];
	char filename2[256];
	char filename3[256];
#ifdef PARTICLES
#ifdef STARFORM
	char filename4[256];
#endif /* STARFORM */
#endif /* PARTICLES */

	int grid_save_flag;
	int particle_save_flag;
	int tracer_save_flag;

	if ( restart_frequency != 0 && step % (2*restart_frequency) == 0 ) {
		grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_BACKUP;
	} else {
		grid_save_flag = particle_save_flag = tracer_save_flag = WRITE_GENERIC;
	}

	if ( particle_output_frequency != 0 && step % particle_output_frequency == 0 ) {
		particle_save_flag = WRITE_SAVE;
	}

	if ( tracer_output_frequency != 0 && step % tracer_output_frequency == 0 ) {
		tracer_save_flag = WRITE_SAVE;
	}

	if ( ( grid_output_frequency != 0 && step % grid_output_frequency == 0 )
#ifdef COSMOLOGY
	     || ( current_output < num_outputs && auni[min_level] >= outputs[current_output] ) 
#endif /* COSMOLOGY */
	     ) {
		/* we're saving information for permenant saving */
		write_restart( WRITE_SAVE, WRITE_SAVE, WRITE_SAVE, NULL );
		current_output++;
	} else if ( restart_frequency != 0 && step % restart_frequency == 0 ) {
		write_restart( grid_save_flag, particle_save_flag, tracer_save_flag, NULL );
	} else {

#ifdef PARTICLES
		if ( particle_save_flag == WRITE_SAVE ) {
			start_time( IO_TIMER );
			start_time( PARTICLE_WRITE_IO_TIMER );

			/* only write out particles, no restart */
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
#ifdef COSMOLOGY
			sprintf( filename1, "%s/%s_a%06.4f.dph", output_directory, jobname, auni[min_level] );
			sprintf( filename2, "%s/%s_a%06.4f.dxv", output_directory, jobname, auni[min_level] );
			sprintf( filename3, "%s/%s_a%06.4f.dpt", output_directory, jobname, auni[min_level] );
#else
			sprintf( filename1, "%s/%s_a%06d.dph", output_directory, jobname, step );
			sprintf( filename2, "%s/%s_a%06d.dxv", output_directory, jobname, step );
			sprintf( filename3, "%s/%s_a%06d.dpt", output_directory, jobname, step );
#endif /* COSMOLOGY */
#else  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
#ifdef COSMOLOGY
			sprintf( filename1, "%s/PMcrda%06.4f.DAT", output_directory, auni[min_level] );
			sprintf( filename2, "%s/PMcrs0a%06.4f.DAT", output_directory, auni[min_level] );
			sprintf( filename3, "%s/pta%06.4f.dat", output_directory, auni[min_level] );
#else
			sprintf( filename1, "%s/PMcrda%06d.DAT", output_directory, step );
			sprintf( filename2, "%s/PMcrs0a%06d.DAT", output_directory, step );
			sprintf( filename3, "%s/pta%06d.dat", output_directory, step );
#endif /* COSMOLOGY */
#endif /* PREFIX_JOBNAME_TO_OUTPUT_FILES */


#ifdef STARFORM
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
#ifdef COSMOLOGY
			sprintf( filename4, "%s/%s_a%06.4f.dst", output_directory, jobname, auni[min_level] );
#else
			sprintf( filename4, "%s/%s_a%06d.dst", output_directory, jobname, step );
#endif /* COSMOLOGY */
#else  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
#ifdef COSMOLOGY
			sprintf( filename4, "%s/stars_a%06.4f.dat", output_directory, auni[min_level] );
#else
			sprintf( filename4, "%s/stars_a%06d.dat", output_directory, step );
#endif /* COSMOLOGY */
#endif /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
			write_particles( filename1, filename2, filename3, filename4 );
#else
			write_particles( filename1, filename2, filename3, NULL );
#endif /* STARFORM */

			end_time( PARTICLE_WRITE_IO_TIMER );
			end_time( IO_TIMER );
		}
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
		if ( tracer_save_flag == WRITE_SAVE ) {
			start_time( IO_TIMER );
			start_time( PARTICLE_WRITE_IO_TIMER );
#ifdef PREFIX_JOBNAME_TO_OUTPUT_FILES
#ifdef COSMOLOGY
			sprintf( filename1, "%s/%s_a%06.4f.dtr", output_directory, jobname, auni[min_level] );
#else
			sprintf( filename1, "%s/%s_a%06d.dtr", output_directory, jobname, step );
#endif /* COSMOLOGY */
#else  /* PREFIX_JOBNAME_TO_OUTPUT_FILES */
#ifdef COSMOLOGY
			sprintf( filename1, "%s/tracers_a%06.4f.dat", output_directory, auni[min_level] );
#else
			sprintf( filename1, "%s/tracers_a%06d.dat", output_directory, step );
#endif /* COSMOLOGY */
#endif /* PREFIX_JOBNAME_TO_OUTPUT_FILES */

			write_hydro_tracers( filename1 );
			end_time( PARTICLE_WRITE_IO_TIMER );
			end_time( IO_TIMER );
		}
#endif /* HYDRO_TRACERS */
	}
}


void set_grid_file_mode(int mode)
{
  if(mode>=-3 && mode<=3) grid_file_mode = mode;
}


void set_particle_file_mode(int mode)
{
  if(mode>=-2 && mode<=2) particle_file_mode = mode;
}


/*
//  Create multiple versions of restart_load_balance using
//  C-style templates
*/
#define FUNCTION                      restart_load_balanceP1
#define PARTICLE_FLOAT                float
#include "io1.def"

#define FUNCTION                      restart_load_balanceP2
#define PARTICLE_FLOAT                double
#include "io1.def"

void restart_load_balance( char *grid_filename, char *particle_header_filename, char *particle_data ) {
  switch(particle_file_mode)
    {
    case 0:
    case 1:
    case -1:
      {
	restart_load_balanceP2(grid_filename,particle_header_filename,particle_data);
	break;
      }
    case 2:
    case -2:
      {
	restart_load_balanceP1(grid_filename,particle_header_filename,particle_data);
	break;
      }
    default:
      {
	cart_error("Invalid particle_file_mode=%d",particle_file_mode);
      }
    }
}

#ifdef PARTICLES

int compare_particle_ids( const void *a, const void *b ) {
	return ( particle_id[*(int *)a] - particle_id[*(int *)b] );
}

void write_particles( char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename ) {
	int i;
	int num_parts;
	char desc[46];
	int processor_heap[MAX_PROCS];
	int *particle_order;
	int *page_ids[MAX_PROCS];
	double *page[MAX_PROCS];
	double *output_page;
	double *times_page[MAX_PROCS];
	double *output_times;
	int count[MAX_PROCS];
	int pos[MAX_PROCS];
	int pages[MAX_PROCS];
	int num_parts_in_page, num_parts_per_page;
	int num_parts_per_proc_page;
	int num_pages, index;
	int local_count;
	int current_id;
	int size;
	FILE *output;
	FILE *timestep_output;
	particle_header header;
	MPI_Status status;
	int left, right, root, smallest;
	int leftproc, rightproc, rootproc, smallestproc;
	int proc;

#ifdef STARFORM
	FILE *stellar_output;
	int num_stars;
	int first_star;
	float *output_stars;
	float *star_page[MAX_PROCS];
#endif /* STARFORM */

	if(particle_file_mode > 0) particle_file_mode = 0; /* Auto-reset */

	num_parts_per_page = num_row*num_row;
	num_parts_per_proc_page = num_parts_per_page/num_procs;
	num_pages = (num_particles_total-1) / num_parts_per_page + 1;

	/* construct mapping of id -> particle index */
	particle_order = cart_alloc(int, num_local_particles );
	num_parts = 0;
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
			particle_order[num_parts++] = i;
		}
	}
	cart_assert( num_parts == num_local_particles );

	qsort( particle_order, num_parts, sizeof(int), compare_particle_ids );

	/* write file header */
	if ( local_proc_id == MASTER_NODE ) {
		output = fopen( header_filename, "w");
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing", header_filename );
		}

		/* set up header and write it out */
#ifdef COSMOLOGY
		header.aunin	= auni[min_level];
		header.auni0	= auni_init;
		header.astep	= abox[min_level] - abox_old[min_level];
		header.h100     = cosmology->h;
		header.OmM0	= cosmology->OmegaM;
		header.OmL0	= cosmology->OmegaL;
		header.OmK0	= cosmology->OmegaK;
		header.OmB0	= cosmology->OmegaB;
		header.DelDC	= cosmology->DeltaDC;
		header.abox	= abox[min_level];
		header.Hbox	= (abox_from_tcode(tl[min_level]+0.5*dtl[min_level])-abox_from_tcode(tl[min_level]-0.5*dtl[min_level]))/dtl[min_level];
#else
		header.aunin	= 1.0;
		header.auni0	= 1.0;
		header.astep	= 0.0;
		header.h100     = 0.0;
		header.OmM0	= primary_units->mass;
		header.OmB0	= primary_units->time;
		header.OmL0	= primary_units->length;
		header.OmK0	= 0.0;
		header.DelDC	= 0.0;
		header.abox	= tl[min_level];
		header.Hbox	= 0.0;
#endif /* COSMOLOGY */
		header.amplt	= 0.0;
		header.istep	= step;
		header.partw	= 0.0;
		header.tintg	= tintg;
		header.ekin		= ekin;
		header.ekin1	= ekin1;
		header.ekin2	= 0.0;
		header.au0		= au0;
		header.aeu0		= aeu0;
		header.Nrow		= num_row;
		header.Ngrid	= num_grid;
		header.Nspecies	= num_particle_species;
		header.Nseed	= 0.0;
		header.Wp5	= 0.0;

		header.magic1	= PARTICLE_HEADER_MAGIC;  /* for indentifying legacy files */
		header.magic2	= PARTICLE_HEADER_MAGIC;  /* for indentifying legacy files */

		for ( i = 0; i < num_particle_species; i++ ) {
			header.mass[i] = particle_species_mass[i];
			header.num[i] = particle_species_indices[i+1];
		}

		/* write jobname to header desc (head of file) */
		snprintf( desc, 45, "%s", jobname );
		for ( i = strlen(jobname); i < 45; i++ ) {
			desc[i] = ' '; /* for fortran version */
		}
		desc[45] = 0;

		size = sizeof(particle_header)+45;
		fwrite( &size, sizeof(int), 1, output );
		fwrite( desc, sizeof(char), 45, output );
		fwrite( &header, sizeof(particle_header), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		fclose( output );

		/* now prepare to write data file */
		output = fopen( data_filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing.", data_filename );
		}

		if ( timestep_filename != NULL ) {
			timestep_output = fopen( timestep_filename, "w" );
			if ( timestep_output == NULL ) {
				cart_error("Unable to open %s for writing.", timestep_filename );
			}

			size = num_particles_total * sizeof(double);
			fwrite( &size, sizeof(int), 1, timestep_output );

			output_times = cart_alloc(double, num_parts_per_page );

			for ( i = 1; i < num_procs; i++ ) {
			        times_page[i] = cart_alloc(double, num_parts_per_proc_page );
			}
		}

		/* allocate space to receive pages of particles from other procs */
		for ( i = 1; i < num_procs; i++ ) {
			page[i] = cart_alloc(double, 2*nDim*num_parts_per_proc_page );
			page_ids[i] = cart_alloc(int, num_parts_per_proc_page );
		}

		/* allocate actual page which will be written */
		output_page = cart_alloc(double, 2*nDim*num_parts_per_page );

		processor_heap[0] = MASTER_NODE;
		page_ids[0] = cart_alloc(int, 1 );
		pos[0] = 0;
		if ( num_local_particles > 0 ) {
			page_ids[0][0] = particle_id[particle_order[0]];
		} else {
			page_ids[0][0] = -1;
		}

		/* receive initial pages */
		for ( i = 1; i < num_procs; i++ ) {
			MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[i] );
			MPI_Recv( page[i], 2*nDim*num_parts_per_proc_page, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
			
			if ( timestep_filename != NULL ) {
				MPI_Recv( times_page[i], num_parts_per_proc_page, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
			}

			pos[i] = 0;
			processor_heap[i] = i;
			pages[i] = 1;

			/* fix for case where one processor has no particles */
			if ( count[i] == 0 ) {
				page_ids[i][0] = -1;
			}
		}

		/* build processor heap ordered by minimum particle id */
		for ( i = num_procs/2-1; i >= 0; i-- ) {
			root = i;
			left = 2*root+1;
			while ( left < num_procs ) {
				leftproc = processor_heap[left];
				rootproc = processor_heap[root];

				if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
						page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
					smallest = left;
				} else {
					smallest = root;
				}

				right = left+1;
				if ( right < num_procs && ( page_ids[processor_heap[smallest]][0] == -1 || 
						( page_ids[processor_heap[right]][0] != -1 && 
						  page_ids[processor_heap[right]][0] < 
						  page_ids[processor_heap[smallest]][0] ) ) ) {
					smallest = right;
				}

				if ( smallest != root ) {
					/* swap */
					processor_heap[root] = processor_heap[smallest];
					processor_heap[smallest] = rootproc;
					root = smallest;
					left = 2*root+1;
				} else {
					left = num_procs;
				}
			}
		}

		current_id = 0;
		local_count = 0;
		for ( i = 0; i < num_pages; i++ ) {
			/* construct page */
			num_parts = 0;

			if ( i == num_pages - 1 ) {
				num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
			} else {
				num_parts_in_page = num_parts_per_page;
			}

			cart_assert( num_parts_in_page > 0 && num_parts_in_page <= num_parts_per_page );

			while ( num_parts < num_parts_in_page ) {
				/* pop minimum processor off the heap */
				proc = processor_heap[0];
				cart_assert( page_ids[proc][pos[proc]] == current_id );

				if ( proc == local_proc_id ) {
					/* add from our list */
					while ( num_parts < num_parts_in_page &&
							local_count < num_local_particles &&
							particle_id[particle_order[local_count]] == current_id ) {
						index = particle_order[local_count];

						output_page[num_parts] = particle_x[index][0] + 1.0;
						output_page[num_parts_per_page+num_parts] = particle_x[index][1] + 1.0;
						output_page[2*num_parts_per_page+num_parts] = particle_x[index][2] + 1.0;
						output_page[3*num_parts_per_page+num_parts] = particle_v[index][0];
						output_page[4*num_parts_per_page+num_parts] = particle_v[index][1];
						output_page[5*num_parts_per_page+num_parts] = particle_v[index][2];

						if ( timestep_filename != NULL ) {
							output_times[num_parts] = particle_dt[index];
						}

						local_count++;
						num_parts++;
						current_id++;
					}

					if ( local_count < num_local_particles ) {
						page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
					} else {
						page_ids[local_proc_id][0] = -1;
					}
				} else {
					while ( num_parts < num_parts_in_page &&
							pos[proc] < count[proc] &&
							page_ids[proc][pos[proc]] == current_id ) {

						output_page[num_parts] = 			page[proc][2*nDim*pos[proc]]   + 1.0;
						output_page[num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+1] + 1.0;
						output_page[2*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+2] + 1.0;
						output_page[3*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+3];
						output_page[4*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+4];
						output_page[5*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+5];

						if ( timestep_filename != NULL ) {
							output_times[num_parts] = times_page[proc][pos[proc]];
						}

						num_parts++;
						current_id++;
						pos[proc]++;

						/* if we've run out, refill page */
						if ( pos[proc] == count[proc] ) { 
							if ( count[proc] == num_parts_per_proc_page ) {
								MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT, 
									proc, pages[proc], MPI_COMM_WORLD, &status );
								MPI_Get_count( &status, MPI_INT, &count[proc] );
								MPI_Recv( page[proc], 2*nDim*num_parts_per_proc_page, 
									MPI_DOUBLE, proc, pages[proc], MPI_COMM_WORLD, &status );

								if ( timestep_filename != NULL ) {
									MPI_Recv( times_page[proc], num_parts_per_proc_page,
											MPI_DOUBLE, proc, pages[proc], MPI_COMM_WORLD, &status );
								}

								pos[proc] = 0;
								pages[proc]++;

								if ( count[proc] == 0 ) {
									page_ids[proc][0] = -1;
								}
							} else {
								pos[proc] = 0;
								count[proc] = 0;
								page_ids[proc][0] = -1;
							}
						}
					}
				}

				/* add back to heap (heapify) */
				root = 0;
				left = 1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][pos[rootproc]] == -1 ||
							( page_ids[leftproc][pos[leftproc]] != -1 && 
							  page_ids[leftproc][pos[leftproc]] < 
							  page_ids[rootproc][pos[rootproc]] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs ) {
						rightproc = processor_heap[right];
						smallestproc = processor_heap[smallest];

						if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
								( page_ids[rightproc][pos[rightproc]] != -1 &&
								  page_ids[rightproc][pos[rightproc]] < 
								  page_ids[smallestproc][pos[smallestproc]] ) ) {
							smallest = right;
						}
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			cart_assert( num_parts == num_parts_in_page );

			/* write page */
			fwrite( output_page, sizeof(double), 2*nDim*num_parts_per_page, output );

			/* write particle timesteps */
			if ( timestep_filename != NULL ) {
				fwrite( output_times, sizeof(double), num_parts_per_page, timestep_output );
			}
		}

		fclose( output );
		cart_free( output_page );

		if ( timestep_filename != NULL ) {
			size = num_particles_total * sizeof(double);
                        fwrite( &size, sizeof(int), 1, timestep_output );

			fclose( timestep_output );
			cart_free( output_times );

			for ( i = 1; i < num_procs; i++ ) {
				cart_free( times_page[i] );
			}
		}

		cart_assert( local_count == num_local_particles );

		for ( i = 1; i < num_procs; i++ ) {
			cart_assert( pos[i] == 0 || pos[i] == count[i] );
			cart_free( page[i] );
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			stellar_output = fopen(stellar_filename, "w");

			if ( stellar_output == NULL ) {
				cart_error("ERROR: unable to open %s for writing.", stellar_filename );
			}

			/* allocate buffers */
			for ( i = 1; i < num_procs; i++ ) {
			        star_page[i] = cart_alloc(float, num_parts_per_proc_page );
			}

			output_stars = cart_alloc(float, num_parts_per_page );

			num_stars = particle_species_num[num_particle_species-1];

			if ( num_stars == 0 ) {
				num_pages = 0;
				first_star = 0;
			} else {
				num_pages = (num_stars-1) / num_parts_per_page + 1;

				first_star = 0;
				while ( first_star < num_local_particles && 
						!particle_is_star( particle_order[first_star] ) ) {
					first_star++;
				}
			}

			/* write header */
			size = 2*sizeof(double);
			fwrite( &size, sizeof(int), 1, stellar_output );
			fwrite( &tl[min_level], sizeof(double), 1, stellar_output );
#ifdef COSMOLOGY
			fwrite( &auni[min_level], sizeof(double), 1, stellar_output );
#else
			fwrite( &dtl[min_level], sizeof(double), 1, stellar_output );
#endif /* COSMOLOGY */
			fwrite( &size, sizeof(int), 1, stellar_output );

			size = sizeof(int);
			fwrite( &size, sizeof(int), 1, stellar_output );
			fwrite( &num_stars, sizeof(int), 1, stellar_output );
			fwrite( &size, sizeof(int), 1, stellar_output );

			size = 2*sizeof(double);
			fwrite( &size, sizeof(int), 1, stellar_output );
			fwrite( &total_stellar_mass, sizeof(double), 1, stellar_output );
			fwrite( &total_stellar_initial_mass, sizeof(double), 1, stellar_output );
			fwrite( &size, sizeof(int), 1, stellar_output );

			/* particle_mass (pw) */
			size = num_stars * sizeof(float);
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
				page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || 
							( page_ids[leftproc][0] != -1 && 
							  page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = particle_mass[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									count[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < 
								  page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

			/* star_initial_mass (pw0) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				processor_heap[i] = i;
				pages[i] = 1;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_initial_mass[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

			/* star_tbirth (tbirth) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
				page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_tbirth[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

#ifdef ENRICH
			/* star_metallicity_II (zstII) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_metallicity_II[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );
#endif /* ENRICH */
#ifdef ENRICH_SNIa
			/* star_metallicity_Ia (zstIa) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_metallicity_Ia[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

#endif /* ENRICH_SNIa */

			cart_free( output_stars );

			for ( i = 1; i < num_procs; i++ ) {
				cart_free( star_page[i] );
			}

			fclose( stellar_output );
		}
#endif /* STARFORM */

		for ( i = 0; i < num_procs; i++ ) {
			cart_free( page_ids[i] );
		}
	} else {
		page[local_proc_id] = cart_alloc(double, 2*nDim*num_parts_per_proc_page );
		page_ids[local_proc_id] = cart_alloc(int, num_parts_per_proc_page );
		pages[local_proc_id] = 0;

		if ( timestep_filename != NULL ) {
		        times_page[local_proc_id] = cart_alloc(double, num_parts_per_proc_page );
		}	

		num_parts = 0;
		do {
			/* pack page */
			local_count = 0;
			for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
				page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
				page[local_proc_id][local_count++] = particle_x[particle_order[num_parts]][0];
				page[local_proc_id][local_count++] = particle_x[particle_order[num_parts]][1];
				page[local_proc_id][local_count++] = particle_x[particle_order[num_parts]][2];
				page[local_proc_id][local_count++] = particle_v[particle_order[num_parts]][0];
				page[local_proc_id][local_count++] = particle_v[particle_order[num_parts]][1];
				page[local_proc_id][local_count++] = particle_v[particle_order[num_parts]][2];

				if ( timestep_filename != NULL ) {
					times_page[local_proc_id][i] = particle_dt[particle_order[num_parts]];
				}

				num_parts++;
			}

			/* send the page */
			MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
			MPI_Send( page[local_proc_id], local_count, MPI_DOUBLE, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );

			if ( timestep_filename != NULL ) {
				MPI_Send( times_page[local_proc_id], i, MPI_DOUBLE, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
			}

			pages[local_proc_id]++;
		} while ( i == num_parts_per_proc_page );

		cart_free( page[local_proc_id] );

		if ( timestep_filename != NULL ) {
			cart_free( times_page[local_proc_id] );
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			first_star = 0;
			while ( first_star < num_local_particles && !particle_is_star( particle_order[first_star] ) ) {
				first_star++;
			}

			star_page[local_proc_id] = cart_alloc(float, num_parts_per_proc_page );

			/* send stellar mass (pw) */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = particle_mass[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );

			/* send star_initial_mass (pw0) */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_initial_mass[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );
	
			/* send star_tbirth */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_tbirth[particle_order[num_parts]];
					num_parts++;
				}
				
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );

#ifdef ENRICH
			/* send metallicities */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_metallicity_II[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );
	
#ifdef ENRICH_SNIa
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_metallicity_Ia[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

			cart_free( star_page[local_proc_id] );
		}
#endif /* STARFORM */

		cart_free( page_ids[local_proc_id] );
	}

	cart_free( particle_order );
}

void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ) {
	int i;
	FILE *input;
	nbody_particle_header nbody_header;
	char desc[46];
	int size;
	char *p;

	*nbody_flag = 0;
	*endian = 0;

	/* read file header */
	input = fopen( header_filename, "r");
	if ( input == NULL ) {
		cart_error("Unable to open particle file %s", header_filename );
	}

	fread( &size, sizeof(int), 1, input );
	fread( desc, sizeof(char), 45, input );
	desc[45] = '\0';

	cart_debug( "Particle header file: %s", desc );

	if ( !control_parameter_is_set("jobname") ) {
		cart_debug("setting jobname to header value");
		strcpy( jobname, desc );

		/* trim spaces from jobname */
		p = jobname + strlen(jobname);
		while (*--p == ' ') *p = '\0';
	}

	if ( size != sizeof(particle_header)+45 ) {
		if ( size == sizeof(nbody_particle_header)+45 ) {
			*nbody_flag = 1;
		} else {
			reorder( (char *)&size, sizeof(int) );

			if ( size != sizeof(particle_header)+45 ) {
				if ( size == sizeof(nbody_particle_header)+45 ) {
					*endian = 1;
					*nbody_flag = 1;
				} else {
					cart_error("Size mismatch in reading particle file header %s (%u vs %u)",
							header_filename, size, sizeof(particle_header)+45 );
				}
			} else {
				*endian = 1;
			}
		}
	}

	if ( *nbody_flag ) {
		cart_debug("USING OLD NBODY FILE FORMAT (hope that's what you wanted...)");

		fread( &nbody_header, sizeof(nbody_particle_header), 1, input );

		header->aunin = nbody_header.aunin;
		header->auni0 = nbody_header.auni0;
		header->amplt = nbody_header.amplt;
		header->astep = nbody_header.astep;
		header->istep = nbody_header.istep;
		header->partw = nbody_header.partw;
		header->tintg = nbody_header.tintg;
		header->ekin = nbody_header.ekin;
		header->ekin1 = nbody_header.ekin1;
		header->ekin2 = nbody_header.ekin2;
		header->au0 = nbody_header.au0;
		header->aeu0 = nbody_header.aeu0;
		header->Nrow = nbody_header.Nrow;
		header->Ngrid = nbody_header.Ngrid;
		header->Nspecies = nbody_header.Nspecies;
		header->Nseed = nbody_header.Nseed;
		header->OmM0 = nbody_header.OmM0;
		header->OmL0 = nbody_header.OmL0;
		header->h100 = nbody_header.h100;
		header->Wp5  = nbody_header.Wp5;
		header->OmK0 = nbody_header.OmK0;

		header->OmB0 = 1.0e-4;

		header->magic1    = nbody_header.magic1;
		header->magic2    = nbody_header.magic2;
		header->DelDC    = nbody_header.DelDC;
		header->abox     = nbody_header.abox;
		header->Hbox     = nbody_header.Hbox;

		for ( i = 0; i < 10; i++ ) {
			header->mass[i] = nbody_header.mass[i];
			header->num[i] = nbody_header.num[i];
		}

	} else {
		fread( header, sizeof(particle_header), 1, input );
	}

	fread( &size, sizeof(int), 1, input );
	fclose(input);

	if ( *endian ) {
		reorder( (char *)&header->aunin, sizeof(float) );
		reorder( (char *)&header->auni0, sizeof(float) );
		reorder( (char *)&header->amplt, sizeof(float) );
		reorder( (char *)&header->astep, sizeof(float) );
		reorder( (char *)&header->istep, sizeof(int) );
		reorder( (char *)&header->partw, sizeof(float) );
		reorder( (char *)&header->tintg, sizeof(float) );
		reorder( (char *)&header->ekin, sizeof(float) );
		reorder( (char *)&header->ekin1, sizeof(float) );
		reorder( (char *)&header->ekin2, sizeof(float) );
		reorder( (char *)&header->au0, sizeof(float) );
		reorder( (char *)&header->aeu0, sizeof(float) );
		reorder( (char *)&header->Nrow, sizeof(int) );
		reorder( (char *)&header->Ngrid, sizeof(int) );
		reorder( (char *)&header->Nspecies, sizeof(int) );
		reorder( (char *)&header->Nseed, sizeof(int) );
		reorder( (char *)&header->OmM0, sizeof(float) );
		reorder( (char *)&header->OmL0, sizeof(float) );
		reorder( (char *)&header->h100, sizeof(float) );
		reorder( (char *)&header->Wp5, sizeof(float) );
		reorder( (char *)&header->OmK0, sizeof(float) );

		if ( !(*nbody_flag) ) {
			reorder( (char *)&header->OmB0, sizeof(float) );
		}

		reorder( (char *)&header->magic1, sizeof(float) );
		reorder( (char *)&header->DelDC, sizeof(float) );
		reorder( (char *)&header->abox, sizeof(float) );
		reorder( (char *)&header->Hbox, sizeof(float) );
		reorder( (char *)&header->magic2, sizeof(float) );

		for ( i = 0; i < header->Nspecies; i++ ) {
			reorder( (char *)&header->mass[i], sizeof(float) );
			reorder( (char *)&header->num[i], sizeof(int) );
		}
	}

	/*
	//  NG: Check for legacy files
	*/
	if(header->magic1!= PARTICLE_HEADER_MAGIC || header->magic2 != PARTICLE_HEADER_MAGIC )
	  {
	    /*
	    //  A legacy file with the garbage in it. Set the DC mode to 0
	    */
	    header->DelDC = 0.0;
	    header->abox = header->aunin;

		cart_debug("Detected legacy format particle file (no DC mode set).");
	  }
}


/*
//  Create multiple versions of read_particles using
//  C-style templates
*/
#define FUNCTION                      read_particlesP1T1
#define PARTICLE_FLOAT                float
#define MPI_PARTICLE_FLOAT            MPI_FLOAT
#define PARTICLE_TIMES_FLOAT          float
#define MPI_PARTICLE_TIMES_FLOAT      MPI_FLOAT
#include "io2.def"

#define FUNCTION                      read_particlesP2T1
#define PARTICLE_FLOAT                double
#define MPI_PARTICLE_FLOAT            MPI_DOUBLE
#define PARTICLE_TIMES_FLOAT          float
#define MPI_PARTICLE_TIMES_FLOAT      MPI_FLOAT
#include "io2.def"

#define FUNCTION                      read_particlesP2T2
#define PARTICLE_FLOAT                double
#define MPI_PARTICLE_FLOAT            MPI_DOUBLE
#define PARTICLE_TIMES_FLOAT          double
#define MPI_PARTICLE_TIMES_FLOAT      MPI_DOUBLE
#include "io2.def"

void read_particles( char *header_filename, char *data_filename, 
			char *timestep_filename, char *stellar_filename,
			int num_sfcs, int *sfc_list ) {
  switch(particle_file_mode)
    {
    case 0:
      {
	read_particlesP2T2(header_filename,data_filename,timestep_filename,stellar_filename,num_sfcs,sfc_list);
	break;
      }
    case 1:
    case -1:
      {
	read_particlesP2T1(header_filename,data_filename,timestep_filename,stellar_filename,num_sfcs,sfc_list);
	break;
      }
    case 2:
    case -2:
      {
	read_particlesP1T1(header_filename,data_filename,timestep_filename,stellar_filename,num_sfcs,sfc_list);
	break;
      }
    default:
      {
	cart_error("Invalid particle_file_mode=%d",particle_file_mode);
      }
    }
}
#endif /* PARTICLES */


#ifdef HYDRO

#ifdef HYDRO_TRACERS
int compare_tracer_ids( const void *a, const void *b ) {
	return ( tracer_id[*(int *)a] - tracer_id[*(int *)b] );
}

void write_hydro_tracers( char *filename ) {
	int i, j, k;
	FILE *output;
	int size;
	int icell, index;
	float adum, ainit;
	float boxh, OmM0, OmL0, OmB0, h100;
	int *tracer_order;
	char label[32];
	int num_tracers_in_page;
	int num_tracers_per_page, num_tracers_per_proc_page, num_pages;
	int num_tracers_written, current_id, local_count;
	double *output_page;
	float *output_vars_page;
	int *page_ids[MAX_PROCS];
	float *vars_page[MAX_PROCS];
	double *page[MAX_PROCS];
	int num_pages_sent;
	int num_pages_received[MAX_PROCS];
	int count[MAX_PROCS];
	int pos[MAX_PROCS];
	MPI_Status status;

	cart_debug("writing hydro tracers");

	num_tracers_per_page = num_tracer_row*num_tracer_row;
	num_tracers_per_proc_page = num_tracers_per_page/num_procs;
	num_pages = (num_tracers_total-1) / num_tracers_per_page + 1;

	/* construct mapping of id -> tracer id */
	tracer_order = cart_alloc(int, num_local_tracers );
	j = 0;
	for ( i = 0; i < num_tracers; i++ ) {
		if ( tracer_id[i] != NULL_TRACER ) {
			tracer_order[j++] = i;
		}
	}
	cart_assert( j == num_local_tracers );

	qsort( tracer_order, num_local_tracers, sizeof(int), compare_tracer_ids );

	if ( local_proc_id == MASTER_NODE ) { 
		output = fopen( filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open hydro tracer file %s.", filename );
		}

		/* write header to file */
		size = 256*sizeof(char);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(jobname, sizeof(char), 256, output );
		fwrite(&size, sizeof(int), 1, output );
		
		/* istep, t, dt, adum, ainit */
#ifdef COSMOLOGY
		adum = auni[min_level];
		ainit = auni_init;
#else
		adum = ainit = 1.0;
#endif

		size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &step, sizeof(int), 1, output );
		fwrite( &tl[min_level], sizeof(double), 1, output );
		fwrite( &dtl[min_level], sizeof(double), 1, output );
#ifdef COSMOLOGY
		fwrite( &adum, sizeof(float), 1, output );
		fwrite( &ainit, sizeof(float), 1, output );
#endif
		fwrite( &size, sizeof(int), 1, output );

#ifdef COSMOLOGY
		/* boxh, Om0, Oml0, Omb0, hubble */
		boxh = box_size;
		OmM0 = cosmology->OmegaM;
		OmL0 = cosmology->OmegaL;
		OmB0 = cosmology->OmegaB;
		h100 = cosmology->h;

		size = 5*sizeof(float);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &boxh, sizeof(float), 1, output );
		fwrite( &OmM0, sizeof(float), 1, output );
		fwrite( &OmL0, sizeof(float), 1, output );
		fwrite( &OmB0, sizeof(float), 1, output );
		fwrite( &h100, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );
#endif /* COSMOLOGY */

		/* num_grid */
		size = sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		size = num_grid;
		fwrite( &size, sizeof(int), 1, output );
		size = sizeof(int);
		fwrite( &size, sizeof(int), 1, output );

		/* number of tracers & page_size */
		size = 2*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &num_tracers_total, sizeof(int), 1, output );
		fwrite( &num_tracer_row, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* hydro variables traced */
		size = sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &num_hydro_vars_traced, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		size = 32*num_hydro_vars_traced*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );

		for ( i = 0; i < num_hydro_vars_traced; i++ ) {
			snprintf( label, 32, "%s", hydro_vars_traced_labels[i] );
			fwrite( label, sizeof(char), 32, output );
		}

		fwrite( &size, sizeof(int), 1, output );

		/* allocate space to receive pages of tracers from other procs */
		for ( i = 1; i < num_procs; i++ ) {
			page_ids[i] = cart_alloc(int, num_tracers_per_proc_page );
			page[i] = cart_alloc(double, nDim*num_tracers_per_proc_page );
			vars_page[i] = cart_alloc(float, num_hydro_vars_traced*num_tracers_per_proc_page );
		}

		/* allocate actual page which will be written */
		output_page = cart_alloc(double, nDim*num_tracers_per_page );

		/* receive initial pages */
		for ( i = 1; i < num_procs; i++ ) {
			num_pages_received[i] = 0;
			MPI_Recv( page_ids[i], num_tracers_per_proc_page, 
				MPI_INT, i, num_pages_received[i], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[i] );
			cart_assert( count[i] >= 0 && count[i] <= num_tracers_per_proc_page );
			MPI_Recv( page[i], nDim*num_tracers_per_proc_page, MPI_DOUBLE, i, 
				num_pages_received[i], MPI_COMM_WORLD, &status );
			pos[i] = 0;
			num_pages_received[i]++;
		}

		current_id = 0;
		local_count = 0;
		for ( i = 0; i < num_pages; i++ ) {
			/* construct page */
			num_tracers_written = 0;
			if ( i == num_pages - 1 ) {
				num_tracers_in_page = num_tracers_total - num_tracers_per_page*(num_pages-1);
			} else {
				num_tracers_in_page = num_tracers_per_page;
			}

			cart_assert( num_tracers_in_page >= 0 && num_tracers_in_page <= num_tracers_per_page );

			while ( num_tracers_written < num_tracers_in_page ) {
				/* add from our list */
				while ( num_tracers_written < num_tracers_in_page &&
						local_count < num_local_tracers &&
						tracer_id[tracer_order[local_count]] == current_id ) {
					index = tracer_order[local_count];
					
					output_page[num_tracers_written]			= tracer_x[index][0] + 1.0;
					output_page[num_tracers_per_page+num_tracers_written] 	= tracer_x[index][1] + 1.0;
					output_page[2*num_tracers_per_page+num_tracers_written]	= tracer_x[index][2] + 1.0;

					local_count++;
					num_tracers_written++;
					current_id++;
				}

				/* add from each processor in turn */
				for ( j = 1; j < num_procs; j++ ) {
					while ( num_tracers_written < num_tracers_in_page &&
							pos[j] < count[j] &&
							page_ids[j][pos[j]] == current_id ) {

						output_page[num_tracers_written]			= page[j][nDim*pos[j]] + 1.0;
						output_page[num_tracers_per_page+num_tracers_written]	= page[j][nDim*pos[j]+1] + 1.0;
						output_page[2*num_tracers_per_page+num_tracers_written] = page[j][nDim*pos[j]+2] + 1.0;

						num_tracers_written++;
						current_id++;
						pos[j]++;

						/* if we've run out, refill page */
						if ( pos[j] == count[j] && count[j] == num_tracers_per_proc_page ) {
							MPI_Recv( page_ids[j], num_tracers_per_proc_page, MPI_INT,
								j, num_pages_received[j], MPI_COMM_WORLD, &status );
							MPI_Get_count( &status, MPI_INT, &count[j] );
							cart_assert( count[j] >= 0 && count[j] <= num_tracers_per_proc_page );
							MPI_Recv( page[j], nDim*num_tracers_per_proc_page,
								MPI_DOUBLE, j, num_pages_received[j], 
								MPI_COMM_WORLD, &status );
							pos[j] = 0;
							num_pages_received[j]++;
						}
					}
				}
			}

			cart_assert( num_tracers_written == num_tracers_in_page );

			/* write page */
			fwrite( output_page, sizeof(double), nDim*num_tracers_per_page, output );
		}

		cart_free( output_page );

                /* allocate actual page which will be written */
		output_vars_page = cart_alloc(float, num_hydro_vars_traced*num_tracers_per_page );

		/* receive initial pages */
		for ( i = 1; i < num_procs; i++ ) {
			MPI_Recv( page_ids[i], num_tracers_per_proc_page, MPI_INT, 
				i, num_pages_received[i], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[i] );
			cart_assert( count[i] >= 0 && count[i] <= num_tracers_per_proc_page );
			MPI_Recv( vars_page[i], num_hydro_vars_traced*num_tracers_per_proc_page, 
				MPI_FLOAT, i, num_pages_received[i], MPI_COMM_WORLD, &status );
			pos[i] = 0;
			num_pages_received[i]++;
		}

		current_id = 0;
		local_count = 0;
		for ( i = 0; i < num_pages; i++ ) {
			/* construct page */
			num_tracers_written = 0;
			if ( i == num_pages - 1 ) {
				num_tracers_in_page = num_tracers_total - num_tracers_per_page*(num_pages-1);
			} else {
				num_tracers_in_page = num_tracers_per_page;
			}

			cart_assert( num_tracers_in_page >= 0 && num_tracers_in_page <= num_tracers_per_page );

			while ( num_tracers_written < num_tracers_in_page ) {
				/* add from our list */
				while ( num_tracers_written < num_tracers_in_page &&
						local_count < num_local_tracers &&
						tracer_id[tracer_order[local_count]] == current_id ) {

					index = tracer_order[local_count];
					icell = cell_find_position( tracer_x[index] );

					for ( j = 0; j < num_hydro_vars_traced; j++ ) {
						cart_assert( j*num_tracers_per_page + num_tracers_written < num_hydro_vars_traced*num_tracers_per_page );
						output_vars_page[j*num_tracers_per_page + num_tracers_written] = 
								cell_vars[icell][ hydro_vars_traced[j] ];
					}

					local_count++;
					num_tracers_written++;
					current_id++;
				}

				/* add from each processor in turn */
				for ( j = 1; j < num_procs; j++ ) {
					while ( num_tracers_written < num_tracers_in_page &&
							pos[j] < count[j] &&
							page_ids[j][pos[j]] == current_id ) {

						for ( k = 0; k < num_hydro_vars_traced; k++ ) {
							output_vars_page[k*num_tracers_per_page + num_tracers_written ] =
								vars_page[j][num_hydro_vars_traced*pos[j]+k];
						}
					
						num_tracers_written++;
						current_id++;
						pos[j]++;

						/* if we've run out, refill page */
						if ( pos[j] == count[j] && count[j] == num_tracers_per_proc_page ) {
							MPI_Recv( page_ids[j], num_tracers_per_proc_page, MPI_INT,
								j, num_pages_received[j], MPI_COMM_WORLD, &status );
							MPI_Get_count( &status, MPI_INT, &count[j] );
							cart_assert( count[j] >= 0 && count[j] <= num_tracers_per_proc_page );
							MPI_Recv( vars_page[j], num_hydro_vars_traced*num_tracers_per_proc_page,
								MPI_FLOAT, j, num_pages_received[j], 
								MPI_COMM_WORLD, &status );
							pos[j] = 0;
							num_pages_received[j]++;
						}
					}
				}
			}

			cart_assert( num_tracers_written == num_tracers_in_page );

			/* write page */
			fwrite( output_vars_page, sizeof(float), num_hydro_vars_traced*num_tracers_per_page, output );
		}		

		cart_free( output_vars_page );

		for ( i = 1; i < num_procs; i++ ) {
			cart_free( vars_page[i] );
			cart_free( page[i] );
			cart_free( page_ids[i] );
		}
		fclose( output );
	} else {
		num_pages_sent = 0;
		page_ids[local_proc_id] = cart_alloc(int, num_tracers_per_proc_page );
		page[local_proc_id] = cart_alloc(double, nDim*num_tracers_per_proc_page );
		vars_page[local_proc_id] = cart_alloc(float, num_hydro_vars_traced*num_tracers_per_proc_page );

		num_tracers_written = 0;
		do {
			/* pack page */
			local_count = 0;
			for ( i = 0; i < num_tracers_per_proc_page && num_tracers_written < num_local_tracers; i++ ) {
				index = tracer_order[num_tracers_written];
				cart_assert( index >= 0 && index < num_tracers );
				cart_assert( tracer_id[index] >= 0 && tracer_id[index] < num_tracers_total );
				cart_assert( i >= 0 && i < num_tracers_per_proc_page );
				page_ids[local_proc_id][i] = tracer_id[index];
				page[local_proc_id][local_count++] = tracer_x[index][0];
				page[local_proc_id][local_count++] = tracer_x[index][1];
				page[local_proc_id][local_count++] = tracer_x[index][2];
				num_tracers_written++;
			}

			cart_assert( i <= num_tracers_per_proc_page );

			/* send the page */
			MPI_Send( page_ids[local_proc_id], i, MPI_INT, 
				MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			MPI_Send( page[local_proc_id], local_count, MPI_DOUBLE, 
				MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			num_pages_sent++;
		} while ( i == num_tracers_per_proc_page );

		cart_free( page[local_proc_id] );

		/* send hydro traced variables */
		num_tracers_written = 0;
		do {
			local_count = 0;
			for ( i = 0; i < num_tracers_per_proc_page && num_tracers_written < num_local_tracers; i++ ) {
				index = tracer_order[num_tracers_written];
				cart_assert( index >= 0 && index < num_tracers );
				cart_assert( tracer_id[index] >= 0 && tracer_id[index] < num_tracers_total );

				page_ids[local_proc_id][i] = tracer_id[index];
				icell = cell_find_position( tracer_x[index] );

				for ( j = 0; j < num_hydro_vars_traced; j++ ) {
					vars_page[local_proc_id][local_count++] = cell_vars[icell][ hydro_vars_traced[j] ];
				}

				num_tracers_written++;
			}

			cart_assert( i <= num_tracers_per_proc_page );

			MPI_Send( page_ids[local_proc_id], i, MPI_INT, 	
				MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			MPI_Send( vars_page[local_proc_id], local_count, 
				MPI_FLOAT, MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			num_pages_sent++;
		} while ( i == num_tracers_per_proc_page );

		cart_free( vars_page[local_proc_id] );
		cart_free( page_ids[local_proc_id] );
	}

	cart_free( tracer_order );
}

void read_hydro_tracers( char *filename ) {
	int i, j, k;
	FILE *input;
	int size, endian;
	char tracerjobname[256];
	double tmpt, tmpdt;
	int tmp_num_grid, tmp_num_hydro_vars_traced;
	int num_read;
	int proc;
	int current_id;
	int tracer;
	int index, icell;
	int coords[nDim];
	float ainit, adum;
	float boxh, Om0, Oml0, Omb0, h;
	int *tracer_order;
	char label[32];
	int num_tracers_in_page;
	int num_tracers_per_page, num_tracers_per_proc_page, num_pages;
	int num_tracers_written, local_count;
	double *input_page;
	double *x, *y, *z;
	int *page_ids[MAX_PROCS];
	double *page[MAX_PROCS];
	int count[MAX_PROCS];
	int current_page[MAX_PROCS];
	int pos[MAX_PROCS];
	MPI_Status status;

	if ( local_proc_id == MASTER_NODE ) {
		input = fopen( filename, "r" );
		if ( input == NULL ) {
			cart_error("Unable to open hydro tracer file %s.", filename );
		}

		endian = 0;

                /* write header to file */
		fread( &size, sizeof(int), 1, input );

		if ( size != 256*sizeof(char) ) {
			reorder( (char *)&size, sizeof(int) );

			if ( size == 256*sizeof(char) ) {
				endian = 1;
				cart_debug("switching endianness of tracer file");
			} else {
				cart_debug("hydro tracer file %s appears to be corrupt!", filename );
			}
		}

		fread( tracerjobname, sizeof(char), 256, input );
		fread(&size, sizeof(int), 1, input );

		cart_debug("tracer job name: %s", tracerjobname );

		fread( &size, sizeof(int), 1, input );
		fread( &step, sizeof(int), 1, input );
		fread( &tmpt, sizeof(double), 1, input );
		fread( &tmpdt, sizeof(double), 1, input );
#ifdef COSMOLOGY
		fread( &adum, sizeof(float), 1, input );
		fread( &ainit, sizeof(float), 1, input );
#endif /* COSMOLOGY */
		fread( &size, sizeof(int), 1, input );

#ifdef COSMOLOGY 
		/* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &Om0, sizeof(float), 1, input );
		fread( &Oml0, sizeof(float), 1, input );
		fread( &Omb0, sizeof(float), 1, input );
		fread( &h, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
#endif /* COSMOLOGY */

		/* num_grid */
		fread( &size, sizeof(int), 1, input );
		fread( &tmp_num_grid, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );	

		if ( endian ) {
			reorder( (char *)&tmp_num_grid, sizeof(int) );
		}

		if ( tmp_num_grid != num_grid ) {
			cart_error("num_grid in %s doesn't match compiled: %u vs %u\n", filename, tmp_num_grid, num_grid );
		}

		/* number of tracers & page_size */
		fread( &size, sizeof(int), 1, input );
		fread( &num_tracers_total, sizeof(int), 1, input );
		fread( &num_tracer_row, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&num_tracers_total, sizeof(int) );
			reorder( (char *)&num_tracer_row, sizeof(int) );
		}

		/* hydro variables traced */
		fread( &size, sizeof(int), 1, input );
		fread( &tmp_num_hydro_vars_traced, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&tmp_num_hydro_vars_traced, sizeof(int) );
		}

		fread( &size, sizeof(int), 1, input );

		for ( i = 0; i < tmp_num_hydro_vars_traced; i++ ) {
			fread( label, sizeof(char), 32, input );
		}

		fread( &size, sizeof(int), 1, input );
	}

	MPI_Bcast( &num_tracer_row, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &num_tracers_total, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	num_tracers_per_page = num_tracer_row*num_tracer_row;
	num_tracers_per_proc_page = num_tracers_per_page/num_procs;
	num_pages = (num_tracers_total-1) / num_tracers_per_page + 1;

	if ( local_proc_id == MASTER_NODE ) {
		input_page = cart_alloc(double, nDim*num_tracers_per_page );

		x = input_page;
		y = &input_page[num_tracers_per_page];
		z = &input_page[2*num_tracers_per_page];

		/* allocate buffer space for tracers on other processors */
		for ( i = 1; i < num_procs; i++ ) {
			page_ids[i] = cart_alloc(int, num_tracers_per_proc_page );
			page[i] = cart_alloc(double, nDim*num_tracers_per_proc_page );
			count[i] = 0;
			current_page[i] = 0;
		}

		current_id = 0;

		for ( i = 0; i < num_pages; i++ ) {
			if ( i == num_pages - 1 ) {
				num_tracers_in_page = num_tracers_total - num_tracers_per_page*(num_pages-1);
			} else {
				num_tracers_in_page = num_tracers_per_page;
			}

			num_read = fread( input_page, sizeof(double), nDim*num_tracers_per_page, input );

			if ( num_read != nDim*num_tracers_per_page ) {
				cart_error("Error reading from tracer file %s: insufficient data", filename );
			}

			if ( endian ) {
				for ( j = 0; j < num_tracers_in_page; j++ ) {
					reorder( (char *)&x[j], sizeof(double) );
					reorder( (char *)&y[j], sizeof(double) );
					reorder( (char *)&z[j], sizeof(double) );
				}
			}
	
			for ( j = 0; j < num_tracers_in_page; j++ ) {
				x[j] -= 1.0;
				y[j] -= 1.0;
				z[j] -= 1.0;

				/* enforce periodic boundary conditions */
				if ( x[j] < 0.0 ) {
					x[j] += (double)num_grid;
				} else if ( x[j] >= (double)num_grid ) {
					x[j] -= (double)num_grid;
				}

				if ( y[j] < 0.0 ) {
					y[j] += (double)num_grid;
				} else if ( y[j] >= (double)num_grid ) {
					y[j] -= (double)num_grid;
				}

				if ( z[j] < 0.0 ) {
					z[j] += (double)num_grid;
				} else if ( z[j] >= (double)num_grid ) {
					z[j] -= (double)num_grid;
				}

				coords[0] = (int)(x[j]);
				coords[1] = (int)(y[j]);
				coords[2] = (int)(z[j]);
	
				index = sfc_index( coords );
				cart_assert( index >= 0 && index < max_sfc_index );
				proc = processor_owner( index );
	
				if ( proc == MASTER_NODE ) {
					tracer = tracer_alloc( current_id );

					tracer_x[tracer][0] = x[j];
					tracer_x[tracer][1] = y[j];
					tracer_x[tracer][2] = z[j];
				} else {
					page_ids[proc][count[proc]] = current_id;
					page[proc][nDim*count[proc]] = x[j];
					page[proc][nDim*count[proc]+1] = y[j];
					page[proc][nDim*count[proc]+2] = z[j];
					count[proc]++;

					if ( count[proc] == num_tracers_per_proc_page ) {
						MPI_Send( page_ids[proc], num_tracers_per_proc_page, MPI_INT, proc,
							current_page[proc], MPI_COMM_WORLD );
						MPI_Send( page[proc], nDim*num_tracers_per_proc_page,
							MPI_DOUBLE, proc, current_page[proc], MPI_COMM_WORLD );
						count[proc] = 0;
						current_page[proc]++;
					}
				}

				current_id++;
			}
		}

		/* send final pages */
		for ( proc = 1; proc < num_procs; proc++ ) {
			MPI_Send( page_ids[proc], count[proc], MPI_INT, proc, current_page[proc], MPI_COMM_WORLD );
			MPI_Send( page[proc], nDim*count[proc], MPI_DOUBLE, proc, current_page[proc], MPI_COMM_WORLD );
			cart_free( page_ids[proc] );
			cart_free( page[proc] );
		}

		cart_free( input_page );

		fclose( input );
	} else {

		page_ids[local_proc_id] = cart_alloc(int, num_tracers_per_proc_page );
		page[local_proc_id] = cart_alloc(double, nDim*num_tracers_per_proc_page );
		count[local_proc_id] = num_tracers_per_proc_page;
		current_page[local_proc_id] = 0;

		while ( count[local_proc_id] == num_tracers_per_proc_page ) {
			MPI_Recv( page_ids[local_proc_id], num_tracers_per_proc_page, MPI_INT,
				MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[local_proc_id] );
			MPI_Recv( page[local_proc_id], nDim*num_tracers_per_proc_page,
				MPI_DOUBLE, MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );

			current_page[local_proc_id]++;

			for ( i = 0; i < count[local_proc_id]; i++ ) {
				tracer = tracer_alloc( page_ids[local_proc_id][i] );

				tracer_x[tracer][0] = page[local_proc_id][nDim*i];
				tracer_x[tracer][1] = page[local_proc_id][nDim*i+1];
				tracer_x[tracer][2] = page[local_proc_id][nDim*i+2];
			}
		}

		cart_free( page[local_proc_id] );
		cart_free( page_ids[local_proc_id] );
	}

	cart_debug("num_local_tracers = %u", num_local_tracers );
	build_tracer_list();
}

#endif /* HYDRO_TRACERS */

void read_gas_ic( char *filename ) {
	int i;
	FILE *input;
	int size;
	float boxh, ainit, astep;
	int ncells;
	int endian;
	int coords[nDim];
	int index;
	int proc, icell;
	MPI_Status status;
	int page_count;
	float *input_page;
	int count[MAX_PROCS];
	float *page[MAX_PROCS];
	int *page_indices[MAX_PROCS];
	int var;
	const int num_gas_vars    = 6;
	const int var_index[] = { 	
		HVAR_GAS_DENSITY, HVAR_MOMENTUM, 
		HVAR_MOMENTUM+1, HVAR_MOMENTUM+2,
		HVAR_GAS_ENERGY, HVAR_INTERNAL_ENERGY };

	if ( local_proc_id == MASTER_NODE ) {
		input = fopen(filename, "r");
		if ( input == NULL ) {
			cart_error("Unable to open %s for reading!", filename );
		}
	
		fread( &size, sizeof(int), 1, input );
		if ( size != sizeof(float) ) {
			reorder( (char *)&size, sizeof(int) );
			endian = 1;
			if ( size != sizeof(float) ) {
				cart_error("Bad file-format in read_cell_ic");
			} 
		} else {
			endian = 0;
		}
	
		fread( &boxh, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		fread( &size, sizeof(int), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &astep, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		fread( &size, sizeof(int), 1, input );
		fread( &ncells, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&ainit, sizeof(float) );
			reorder( (char *)&astep, sizeof(float) );
			reorder( (char *)&ncells, sizeof(int) );
		}

		box_size = boxh;
		MPI_Bcast( &box_size, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

#ifdef COSMOLOGY
		auni_init = ainit;
		MPI_Bcast( &auni_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
#endif /* COSMOLOGY */
	
		cart_debug("boxh = %f", boxh );
		cart_debug("ainit = %f", ainit );
		cart_debug("astep = %f", astep );
		cart_debug("ncells = %u", ncells );

		if ( ncells != num_root_cells ) {
			cart_error("ncells in %s does not match num_root_cells (%u vs %u)", 
					filename, ncells, num_root_cells );
		}

		input_page = cart_alloc(float, num_grid*num_grid );
                for ( proc = 1; proc < num_procs; proc++ ) {
			page[proc] = cart_alloc(float, num_grid*num_grid );
			page_indices[proc] = cart_alloc(int, num_grid*num_grid );
                }

		for ( var = 0; var < num_gas_vars; var++ ) {
			for ( proc = 1; proc < num_procs; proc++ ) {
				count[proc] = 0;
			}

			fread( &size, sizeof(int), 1, input );
			for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				size = fread( input_page, sizeof(float), num_grid*num_grid, input );
				if ( size != num_grid*num_grid ) {
					cart_error("Error reading from file %s", filename );
				}

				if ( endian ) {
					for ( i = 0; i < num_grid*num_grid; i++ ) {
						reorder( (char *)&input_page[i], sizeof(float) );
					}
				}
	
				page_count = 0;
				for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
					for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
						index = sfc_index( coords );
						proc = processor_owner( index );
	
						if ( proc == local_proc_id ) {
							icell = root_cell_location( index );
							cell_var( icell, var_index[var] ) = input_page[page_count];
						} else {
							/* add cell to send page */
							page[proc][count[proc]] = input_page[page_count];
							page_indices[proc][count[proc]] = index;
							count[proc]++;
	
							if ( count[proc] == num_grid*num_grid ) {
								MPI_Send( page[proc], num_grid*num_grid, MPI_FLOAT, 
									proc, 0, MPI_COMM_WORLD );
								MPI_Send( page_indices[proc], num_grid*num_grid, 
									MPI_INT, proc, 0, MPI_COMM_WORLD );
								count[proc] = 0;
							}
						}
	
						page_count++;
					}
				}
			}
			fread( &size, sizeof(int), 1, input );
			
			/* send last variables */
			for ( proc = 1; proc < num_procs; proc++ ) {
				MPI_Send( page[proc], count[proc], MPI_FLOAT, proc, 0, MPI_COMM_WORLD );
				MPI_Send( page_indices[proc], count[proc], MPI_INT, proc, 0, MPI_COMM_WORLD );
			}
		}
	
        	fclose(input);

		cart_free( input_page );
		for ( proc = 1; proc < num_procs; proc++ ) {
			cart_free( page[proc] );
			cart_free( page_indices[proc] );
		}
	} else {
		MPI_Bcast( &box_size, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
#ifdef COSMOLOGY
		MPI_Bcast( &auni_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
#endif /* COSMOLOGY */

		page[local_proc_id] = cart_alloc(float, num_grid*num_grid );
		page_indices[local_proc_id] = cart_alloc(int, num_grid*num_grid );

		for ( var = 0; var < num_gas_vars; var++ ) {
			page_count = num_grid*num_grid;
			while ( page_count == num_grid*num_grid ) {
				MPI_Recv( page[local_proc_id], num_grid*num_grid, MPI_FLOAT, MASTER_NODE, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_FLOAT, &page_count );
				MPI_Recv( page_indices[local_proc_id], num_grid*num_grid, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD, &status );

				for ( i = 0; i < page_count; i++ ) {
					icell = root_cell_location( page_indices[local_proc_id][i] );
					cell_var( icell, var_index[var] ) = page[local_proc_id][i];
				}
			}
		}

		cart_free( page[local_proc_id] );
		cart_free( page_indices[local_proc_id] );
	}

	/* set gas gamma on root level */
	for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
		cell_gas_gamma(i) = constants->gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		cell_electron_internal_energy(i) = cell_gas_internal_energy(i)*constants->wmu/constants->wmu_e;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
	}

	for ( i = 0; i < nDim; i++ ) {
		refinement_volume_min[i] = 0.0;
		refinement_volume_max[i] = num_grid;
	}
}

#endif /* HYDRO */


/*
// NG: Slightly re-written form of grid binary I/O that (a) eliminates 
// code duplication and (b) allows complete flexibility in what is written
// to disk (the latter is needed for adding RT block I/O).
*/


/* two helpers */
void write_grid_binary_top_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells);
void write_grid_binary_lower_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells, int level, int current_level_count, int *current_level);


void write_grid_binary( char *filename ) {
	int i, j, k;
	int size;
	FILE *output;
	float adum, ainit;
	float boxh, OmM0, OmL0, OmB0, h100;
	int minlevel, maxlevel;
	int nextras;
	int *cellrefined;
	int *order;
	int *current_level;
	int proc;
	int level;
	int icell, ioct;
	int page_size;
	int ncell0, sfc_order;
	int current_level_count;
	int next_level_count;
	int page;
	int page_count;
	long total_cells[max_level-min_level+1];
	MPI_Status status;
	char parallel_filename[256];
	int proc_num_cells[MAX_PROCS*(max_level-min_level+1)];
	int file_index, file_parent, file_num_procs;

	int hydro_vars[num_hydro_vars+1];
	int num_other_vars = 0;
	int *other_vars = 0;

	if(grid_file_mode > 0) grid_file_mode = 0; /* Auto-reset */

#ifdef HYDRO
	/*
	// Maintain the same order as in a previous version
	*/
	hydro_vars[0] = HVAR_GAS_DENSITY;
	hydro_vars[1] = HVAR_GAS_ENERGY;
	hydro_vars[2] = HVAR_MOMENTUM + 0;
	hydro_vars[3] = HVAR_MOMENTUM + 1;
	hydro_vars[4] = HVAR_MOMENTUM + 2;
	hydro_vars[5] = HVAR_PRESSURE;
	hydro_vars[6] = HVAR_GAMMA;
	hydro_vars[7] = HVAR_INTERNAL_ENERGY;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	hydro_vars[8] = HVAR_ELECTRON_INTERNAL_ENERGY;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for(j=0; j<num_chem_species; j++)
	  {
	    hydro_vars[num_hydro_vars-num_chem_species+j] = HVAR_ADVECTED_VARIABLES+j;
	  }
#endif /* HYDRO */

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
#ifdef GRAVITY
	num_other_vars++;
#ifdef HYDRO
	num_other_vars++;
#endif /* HYDRO */
#endif
#ifdef RADIATIVE_TRANSFER
	num_other_vars += rt_num_disk_vars;
#endif

	other_vars = cart_alloc(int, num_other_vars );

#ifdef GRAVITY
	other_vars[0] = VAR_POTENTIAL;
	k = 1;
#ifdef HYDRO
	other_vars[1] = VAR_POTENTIAL_HYDRO;
	k = 2;
#endif /* HYDRO */
#else
	k = 0;
#endif
#ifdef RADIATIVE_TRANSFER
	for(j=0; j<rt_num_disk_vars; j++) other_vars[k+j] = rt_disk_offset + j; 
#endif
#endif /* GRAVITY || RADIATIVE_TRANSFER */


	/* ensure consistency of num_output_files */
	num_output_files = min( num_output_files, num_procs );
	num_output_files = max( num_output_files, 1 );

	page_size = num_grid*num_grid;

	/* determine parallel output options */
	file_index = local_proc_id * num_output_files / num_procs;

	file_parent = local_proc_id;
	while ( file_parent > 0 && 
			(file_parent-1)*num_output_files / num_procs == file_index ) {
		file_parent--;
	}	

	file_num_procs = 1;
	while ( file_parent + file_num_procs < num_procs &&
			(file_parent+file_num_procs)*num_output_files / num_procs == file_index ) {
		file_num_procs++;
	}	

	cellrefined = cart_alloc(int, page_size );

	/* get number of cells on each level */
	MPI_Allgather( num_cells_per_level, max_level-min_level+1, MPI_INT,
		proc_num_cells, max_level-min_level+1, MPI_INT, MPI_COMM_WORLD );

	minlevel = min_level;
	maxlevel = max_level_now_global(MPI_COMM_WORLD);

	/* open file handle if parent of parallel file */
	if ( local_proc_id == file_parent ) {
		if ( num_output_files == 1 ) {
			output = fopen(filename,"w");
		} else {
			sprintf( parallel_filename, "%s.%03u", filename, file_index );
			output = fopen(parallel_filename, "w");
		}

		if ( output == NULL ) {
			cart_error( "Unable to open file %s for writing!", filename );
		}

		for ( level = min_level; level <= max_level; level++ ) {
			total_cells[level] = 0;
			for ( proc = local_proc_id; proc < local_proc_id+file_num_procs; proc++ ) {
				total_cells[level] += proc_num_cells[(max_level-min_level+1)*proc+level];
			}
		}
	}

	/* only write one copy of the header information */
	if ( local_proc_id == MASTER_NODE ) {
		size = 256*sizeof(char);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(jobname, sizeof(char), 256, output );
		fwrite(&size, sizeof(int), 1, output );

		/* istep, t, dt, adum, ainit */
#ifdef COSMOLOGY
		adum = auni[min_level];
		ainit = auni_init;
#else
		adum = ainit = 1.0;
#endif /* COSMOLOGY */
		size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &step, sizeof(int), 1, output );
		fwrite( &tl[min_level], sizeof(double), 1, output );
		fwrite( &dtl[min_level], sizeof(double), 1, output );
		fwrite( &adum, sizeof(float), 1, output );
		fwrite( &ainit, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* boxh, OmM0, OmL0, OmB0, hubble */
		boxh = box_size;
#ifdef COSMOLOGY
		OmM0 = cosmology->OmegaM;
		OmL0 = cosmology->OmegaL;
		OmB0 = cosmology->OmegaB;
		h100 = cosmology->h;
#else
		OmM0 = primary_units->mass;
		OmB0 = primary_units->time;
		OmL0 = primary_units->length;
		h100 = 0.0; /* indicator that we have units rather than cosmological parameters */
#endif /* COSMOLOGY */
		size = 5*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &boxh, sizeof(float), 1, output );
		fwrite( &OmM0, sizeof(float), 1, output );
		fwrite( &OmL0, sizeof(float), 1, output );
		fwrite( &OmB0, sizeof(float), 1, output );
		fwrite( &h100, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* nextra (no evidence extras are used...) extra lextra */
		size = sizeof(int);
		nextras = 0;

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &nextras, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* extra */
		size = nextras * sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* lextra */
		size = nextras * 256 * sizeof(char);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );
	
		/* Minlevel, MaxLevelNow */
		size = 2 * sizeof(int);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(&minlevel, sizeof(int), 1, output );
		fwrite(&maxlevel, sizeof(int), 1, output );
		fwrite(&size, sizeof(int), 1, output );

		size = (maxlevel-minlevel+1) * sizeof(double);

		/* tl */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( tl, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* dtl */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( dtl, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* tl_old */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( tl_old, sizeof(double), maxlevel-minlevel+1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* dtl_old */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( dtl_old, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* iSO */
		size = (maxlevel-minlevel+1) * sizeof(int);

#ifdef HYDRO
		fwrite( &size, sizeof(int), 1, output );
		fwrite( level_sweep_dir, sizeof(int), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );
#endif /* HYDRO */

		/* sfc ordering used */
		sfc_order = SFC;
		size = sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &sfc_order, sizeof(int), 1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* refinement volume */
		size = 2*nDim*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( refinement_volume_min, sizeof(float), nDim, output );
		fwrite( refinement_volume_max, sizeof(float), nDim, output );
		fwrite( &size, sizeof(int), 1, output );

#ifdef STARFORM
		/* refinement volume */
		size = 2*nDim*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( star_formation_volume_min, sizeof(float), nDim, output );
		fwrite( star_formation_volume_max, sizeof(float), nDim, output );
		fwrite( &size, sizeof(int), 1, output );
#endif /* STARFORM */

		/* ncell0 */
		ncell0 = num_grid*num_grid*num_grid;
		size = sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &ncell0, sizeof(int), 1, output);
		fwrite( &size, sizeof(int), 1, output );
	}

	/* now start writing pages of root level cell children */
	if ( local_proc_id == file_parent ) {
		size = total_cells[min_level] * sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
	}

	/* holds list of next level octs to write */
	order = cart_alloc(int, (num_cells_per_level[min_level+1] / num_children) );
	next_level_count = 0;
	current_level_count = 0;
	page = 0;

	while ( current_level_count < num_cells_per_level[min_level] ) {
		page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );

		for ( i = 0; i < page_count; i++ ) {
			if ( cell_is_refined(current_level_count) ) {
				cellrefined[i] = tree_cell_count(current_level_count);
				order[next_level_count++] = cell_child_oct[current_level_count];
			} else {
				cellrefined[i] = 1;
			}
			current_level_count++;
		}

		if ( local_proc_id == file_parent ) {
			fwrite( cellrefined, sizeof(int), page_count, output );
		} else {
			MPI_Send( cellrefined, page_count, MPI_INT, file_parent, page, MPI_COMM_WORLD );
			page++;
		}
	}

	if ( local_proc_id == file_parent ) {
		for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
			i = 0;
			page = 0;
			while ( i < proc_num_cells[(max_level-min_level+1)*proc+min_level] ) {
				page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );
				MPI_Recv( cellrefined, page_count, MPI_INT, proc, page, MPI_COMM_WORLD, &status );
				fwrite( cellrefined, sizeof(int), page_count, output );
				i += page_count;
				page++;
			}
		}

		fwrite( &size, sizeof(int), 1, output );
	}

	/* now write pages of root level hydro variables */
	write_grid_binary_top_level_vars(num_hydro_vars,hydro_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells);

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
	write_grid_binary_top_level_vars(num_other_vars,other_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells);
#endif /* GRAVITY || RADIATIVE_TRANSFER */

	/* then write each level's cells in turn */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		if ( local_proc_id == file_parent ) {
			/* write size */
			size = sizeof(long);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &total_cells[level], sizeof(long), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			size = total_cells[level] * sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
		}

		current_level_count = next_level_count;
		cart_assert( num_cells_per_level[level] == current_level_count*num_children );
		next_level_count = 0;
		current_level = order;

		if ( level < maxlevel ) {
			order = cart_alloc(int, ( num_cells_per_level[level+1]/num_children) );
		} else {
			order = cart_alloc(int, 0 );
		}

		i = 0;
		page = 0;
		while ( i < current_level_count ) {
			page_count = min( page_size, num_cells_per_level[level] - i*num_children );

			j = 0;
			while ( j < page_count) {
				ioct = current_level[i];

				for ( k = 0; k < num_children; k++ ) {
					icell = oct_child( ioct, k );

					if ( cell_is_refined(icell) ) {
						cellrefined[j++] = 1;
						order[next_level_count++] = cell_child_oct[icell];
					} else {
						cellrefined[j++] = 0;
					}
				}

				i++;
			}

			if ( local_proc_id == file_parent ) {
				fwrite( cellrefined, sizeof(int), page_count, output );
			} else {
				MPI_Send( cellrefined, page_count, MPI_INT, file_parent, page, MPI_COMM_WORLD );
				page++;
			}
		}

		if ( local_proc_id == file_parent ) {
			for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
				i = 0;
				page = 0;
				while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] ) {
					page_count = min( page_size, 
							proc_num_cells[(max_level-min_level+1)*proc+level] - i );
					MPI_Recv( cellrefined, page_count, MPI_INT, proc, page, MPI_COMM_WORLD, &status );
					fwrite( cellrefined, sizeof(int), page_count, output );
					i += page_count;
					page++;
				}
			}

			fwrite( &size, sizeof(int), 1, output );
		}

		/* now write pages of lower level hydro variables */
		write_grid_binary_lower_level_vars(num_hydro_vars,hydro_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells,level,current_level_count,current_level);

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		write_grid_binary_lower_level_vars(num_other_vars,other_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells,level,current_level_count,current_level);
#endif /* GRAVITY || RADIATIVE_TRANSFER */

		cart_free( current_level );
	}

	if ( local_proc_id == file_parent ) {
		fclose( output );
	}

	if(other_vars != 0) cart_free(other_vars);

	cart_free( order );
	cart_free( cellrefined );

#ifdef RADIATIVE_TRANSFER
	/* Save RF data */
	rtWriteRadiationFieldData(filename,1);
#endif
}


void write_grid_binary_top_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells)
{
  int i, j;
  int size;
  float *cellvars;
  int current_level_count;
  int page, page_count;
  int icell;
  int proc;
  MPI_Status status;

  if(num_out_vars < 1) return;

  cellvars = cart_alloc(float, num_out_vars * page_size );

  /* now write pages of root level hydro variables */
  if ( local_proc_id == file_parent )
    {
      size = num_out_vars * total_cells[min_level] * sizeof(float);
      fwrite( &size, sizeof(int), 1, output );
    }

  current_level_count = 0;
  page = 0;
  while ( current_level_count < num_cells_per_level[min_level] )
    {
      page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );
    
      i = 0;
      while ( i < num_out_vars*page_count )
	{
	  icell = current_level_count;

	  for(j=0; j<num_out_vars; j++)
	    {
	      cellvars[i++] = cell_var(icell,out_var[j]);
	    }

	  current_level_count++;
	}

      if ( local_proc_id == file_parent )
	{
	  fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	} 
      else
	{
	  MPI_Send( cellvars, num_out_vars*page_count, MPI_FLOAT, file_parent, page, MPI_COMM_WORLD );
	  page++;
	}
    }

  if ( local_proc_id == file_parent )
    {
      for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ )
	{
	  i = 0;
	  page = 0;
	  while ( i <  proc_num_cells[(max_level-min_level+1)*proc+min_level] )
	    {
	      page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );
	
	      MPI_Recv( cellvars, num_out_vars*page_count, MPI_FLOAT, proc, page, MPI_COMM_WORLD, &status );
	      fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	      i += page_count;
	      page++;
	    }
	}
      
      fwrite( &size, sizeof(int), 1, output );
    }

  cart_free( cellvars );
}


void write_grid_binary_lower_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells, int level, int current_level_count, int *current_level)
{
  int i, j, k, m;
  int size;
  float *cellvars;
  int page, page_count;
  int icell, ioct;
  int proc;
  MPI_Status status;

  if(num_out_vars < 1) return;

  cellvars = cart_alloc(float, num_out_vars * page_size );

  /* now write pages of root level hydro variables */
  if ( local_proc_id == file_parent )
    {
      size = num_out_vars * total_cells[level] * sizeof(float);
      fwrite( &size, sizeof(int), 1, output );
    }

  i = 0;
  page = 0;
  while ( i < current_level_count )
    {
      page_count = min( page_size, num_cells_per_level[level] - i*num_children );
      
      j = 0;
      while ( j < num_out_vars*page_count ) { 
	ioct = current_level[i];

	for ( k = 0; k < num_children; k++ ) {
	  icell = oct_child( ioct, k );

	  for ( m = 0; m < num_out_vars; m++ )
	    {
	      cellvars[j++] = cell_var(icell,out_var[m]);
	    }
	  
	}
	i++;
      }

      if ( local_proc_id == file_parent )
	{
	  fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	}
      else
	{
	  MPI_Send( cellvars, num_out_vars*page_count, MPI_FLOAT, file_parent, page, MPI_COMM_WORLD );
	  page++;
	}
    }

  if ( local_proc_id == file_parent )
    {
      for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ )
	{
	  i = 0;
	  page = 0;
	  while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] )
	    {
	      page_count = min( page_size, proc_num_cells[(max_level-min_level+1)*proc+level] - i );
	      MPI_Recv( cellvars, num_out_vars*page_count, MPI_FLOAT, proc, page, MPI_COMM_WORLD, &status );
	      fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	      i += page_count;
	      page++;
	    }
	}
      
      fwrite( &size, sizeof(int), 1, output );
    }

  cart_free( cellvars );
}


/* two helpers */
void read_grid_binary_top_level_vars(int num_in_vars, int jskip, int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, int local_file_root_cells, int page_size, int *proc_num_cells, long *proc_cell_index, int *file_sfc_index);
void read_grid_binary_lower_level_vars(int num_in_vars, int jskip, int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, long *total_cells, int page_size, int *proc_num_cells, int level, long *first_page_count, long *proc_first_index, long *proc_cell_index, int *current_level);


void read_grid_binary( char *filename ) {
	int i, j;
	int size;
	int num_read, flag;
	FILE *input;
	char job[256];
	int minlevel, maxlevel;
	double t, dt;
	float adum, ainit;
	float boxh, OmM0, OmL0, OmB0, h100;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int *cellrefined[MAX_PROCS], *cellrefinedbuffer;
	int *order;
	int *current_level;
	int endian;
	int proc;
	int ret;
	int level;
	int icell, ioct;
	int cell_counts;
	int page_size;
	int ncell0, sfc_order;
	long current_level_count, current_read_count;
	long next_level_count;
	int file_index, file_parent;
	int local_file_root_cells;
	int page_count;
	long count, start;
	long total_cells[max_level-min_level+1];
	int file_parent_proc[MAX_PROCS];
	int proc_num_cells[MAX_PROCS];
	int proc_cur_cells[MAX_PROCS];
	int proc_page_count[MAX_PROCS];
	long first_page_count[MAX_PROCS];
	long proc_next_level_octs[MAX_PROCS];
	long proc_cell_index[MAX_PROCS];
	long proc_first_index[MAX_PROCS+1];
	long next_proc_index[MAX_PROCS+1];
	int file_sfc_index[MAX_PROCS+1];
	char parallel_filename[256];
	int num_requests, num_send_requests;
	int continue_reading, ready_to_read;
	MPI_Request requests[2*MAX_PROCS];
	MPI_Request send_requests[MAX_PROCS];
	MPI_Status status;

#ifdef COSMOLOGY
	struct CosmologyParameters temp_cosmo;
#endif /* COSMOLOGY */

	int skip_hvar = -1;
	int skip_ovar = -1;
	int hydro_vars[num_hydro_vars+1];
	int num_other_vars = 0;
	int *other_vars = 0;
	int num_in_hydro_vars, num_in_other_vars;
#ifdef RADIATIVE_TRANSFER
	int rt_var_offset = 0;
#endif

#ifdef HYDRO
	/*
	// Maintain the same order as in a previous version
	 */
	hydro_vars[0] = HVAR_GAS_DENSITY;
	hydro_vars[1] = HVAR_GAS_ENERGY;
	hydro_vars[2] = HVAR_MOMENTUM + 0;
	hydro_vars[3] = HVAR_MOMENTUM + 1;
	hydro_vars[4] = HVAR_MOMENTUM + 2;
	hydro_vars[5] = HVAR_PRESSURE;
	hydro_vars[6] = HVAR_GAMMA;
	hydro_vars[7] = HVAR_INTERNAL_ENERGY;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	hydro_vars[8] = HVAR_ELECTRON_INTERNAL_ENERGY;
#endif
	for(j=0; j<num_chem_species; j++) {
        hydro_vars[num_hydro_vars-num_chem_species+j] = HVAR_ADVECTED_VARIABLES+j;
	}
#endif /* HYDRO */

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
#ifdef GRAVITY
	num_other_vars++;
#ifdef HYDRO
	num_other_vars++;
#endif /* HYDRO */
#endif
#ifdef RADIATIVE_TRANSFER
	rt_var_offset = num_other_vars;
	num_other_vars += rt_num_disk_vars;
#endif

	other_vars = cart_alloc(int, num_other_vars );

#ifdef GRAVITY
	other_vars[0] = VAR_POTENTIAL;
#ifdef HYDRO
	other_vars[1] = VAR_POTENTIAL_HYDRO;
#endif /* HYDRO */
#endif
#ifdef RADIATIVE_TRANSFER
	for(j=0; j<rt_num_disk_vars; j++) other_vars[rt_var_offset+j] = rt_disk_offset + j; 
#endif
#endif /* GRAVITY || RADIATIVE_TRANSFER */

	switch(grid_file_mode)
	  {
	  case  0:
	    {
	      num_in_hydro_vars = num_hydro_vars;
	      num_in_other_vars = num_other_vars;
	      break;
	    }
	  case  1:
	  case -1:
	    {
	      num_in_hydro_vars = num_hydro_vars + 6;
	      num_in_other_vars = num_other_vars + 6;
	      break;
	    }
	  case  2:
	  case -2:
	    {
	      num_in_hydro_vars = num_hydro_vars + 6;
	      num_in_other_vars = num_other_vars + 8;
	      break;
	    }
#ifdef BLASTWAVE_FEEDBACK
	  case  3:
	  case -3:
	    {
	      num_in_hydro_vars = num_hydro_vars - 1;
	      num_in_other_vars = num_other_vars ;
	      skip_hvar = HVAR_BLASTWAVE_TIME - HVAR_GAS_DENSITY; //fundhydrovars
	      break;
	    }
#endif /* BLASTWAVE_FEEDBACK */
	  default:
	    {
	      cart_error("Invalid grid_file_mode=%d in a call to  read_grid_binary",grid_file_mode);
	    }
	  }

	cart_assert( num_output_files >= 1 && num_output_files <= num_procs );

	page_size = num_grid*num_grid;

	/* set up global file information */
	proc = 0;	
	for ( i = 0; i < num_output_files; i++ ) {
		while ( proc*num_output_files / num_procs != i ) {
			proc++;
		}
		file_parent_proc[i] = proc;
	}

	/* determine parallel output options */
	file_index = local_proc_id * num_output_files / num_procs;
	file_parent = file_parent_proc[file_index];

	/* open file handle if parent of parallel file */
	if ( local_proc_id == file_parent ) {
		if ( num_output_files == 1 ) {
			input = fopen(filename,"r");
		} else {
			sprintf( parallel_filename, "%s.%03u", filename, file_index );
			input = fopen(parallel_filename, "r");
		}

		if ( input == NULL ) {
			cart_error( "read_grid_binary: unable to open file %s for reading!", filename );
		}
	}

	/* the header exists only in the first file */
	if ( local_proc_id == MASTER_NODE ) {
		cart_assert( local_proc_id == file_parent );

		fread(&size, sizeof(int), 1, input );
		endian = 0;
		if ( size != 256 ) {
			reorder( (char *)&size, sizeof(int) );
			if ( size != 256 ) {
				cart_error("Error: file %s is corrupted", filename );
			} else {
				endian = 1;
				cart_debug("Reordering bytes (file endianness is opposite program)");
			}
		}

		fread(job, sizeof(char), 256, input );
		if ( !control_parameter_is_set("jobname") ) {
		  cart_debug("setting jobname to header value");
		  strcpy( jobname, job );
		}
		fread(&size, sizeof(int), 1, input );

		/* istep, t, dt, adum, ainit */
		fread( &size, sizeof(int), 1, input );
		fread( &step, sizeof(int), 1, input );
		fread( &t, sizeof(double), 1, input );
		fread( &dt, sizeof(double), 1, input );
		fread( &adum, sizeof(float), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
		
		if ( endian ) {
			reorder( (char *)&step, sizeof(int) );
			reorder( (char *)&ainit, sizeof(float) );
		}

		/* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &OmM0, sizeof(float), 1, input );
		fread( &OmL0, sizeof(float), 1, input );
		fread( &OmB0, sizeof(float), 1, input );
		fread( &h100, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&OmM0, sizeof(float) );
			reorder( (char *)&OmL0, sizeof(float) );
			reorder( (char *)&OmB0, sizeof(float) );
			reorder( (char *)&h100, sizeof(float) );
		}

		box_size = boxh;

#ifdef COSMOLOGY
		auni_init = ainit;

		cosmology_set(OmegaM,OmM0);
		cosmology_set(OmegaL,OmL0);
		cosmology_set(OmegaB,OmB0);
		cosmology_set(h,h100);

		/*
		//  NG: we do not set the DC mode here since we
		//  assume it will be set by the particle reader.
		//  The DC mode is only relevant for cosmology, and 
		//  it is unlikely that a cosmology run will not 
		//  include particles. And if it even does, then 
		//  there is no sense whatsoever to set 
		//  a non-trivial DC mode.
		//
		//  cosmology_set(DeltaDC,0.0);
		 */

		temp_cosmo = *cosmology;

#else
		/* do nothing, units now set for all procs together */ 
#endif /* COSMOLOGY */

		/* nextra (no evidence extras are used...) extra lextra */
		fread( &size, sizeof(int), 1, input );
		fread( &nextras, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&nextras, sizeof(int) );
		}

		/* extra */
		fread( &size, sizeof(int), 1, input );
		fread( extra, sizeof(float), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* lextra */
		fread( &size, sizeof(int), 1, input );
		fread( lextra, 256*sizeof(char), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* Minlevel, MaxLevelNow */
		fread( &size, sizeof(int), 1, input );
		fread( &minlevel, sizeof(int), 1, input );
		fread( &maxlevel, sizeof(int), 1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&minlevel, sizeof(int) );
			reorder( (char *)&maxlevel, sizeof(int) );
		}

		if ( maxlevel > max_level ) {
			cart_error("File %s has more levels than compiled program (%u)", filename, maxlevel );
		}

		cart_assert( minlevel == min_level );

		/* tl */
		fread( &size, sizeof(int), 1, input );
		fread( tl, sizeof(double), maxlevel-minlevel+1, input );
		fread( &size, sizeof(int), 1, input);

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl[i], sizeof(double) );
			}
		}

		/* dtl */
		fread( &size, sizeof(int), 1, input );
		fread( dtl, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl[i], sizeof(double) );
			}
		}

		/* tl_old */
		fread( &size, sizeof(int), 1, input );
		fread( tl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl_old[i], sizeof(double) );
			}
		}

		/* dtl_old */
		fread( &size, sizeof(int), 1, input );
		fread( dtl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );	

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl_old[i], sizeof(double) );
			}
		}

		/* iSO (sweep direction for flux solver) */
#ifdef HYDRO
		fread( &size, sizeof(int), 1, input );
		fread( level_sweep_dir, sizeof(int), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&level_sweep_dir[i], sizeof(int) );
			}
		}

		for ( i = maxlevel+1; i <= max_level; i++ ) {
			level_sweep_dir[i] = 0;
		}
#endif /* HYDRO */

		/* sfc ordering used */
		fread( &size, sizeof(int), 1, input );
		fread( &sfc_order, sizeof(int), 1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&sfc_order, sizeof(int) );
		}

		if ( sfc_order != SFC ) {
			cart_error("File has different sfc indexing (%d) than program (%d)",sfc_order,SFC);
		}

		/* refinement volume */
		fread( &size, sizeof(int), 1, input );
		fread( refinement_volume_min, sizeof(float), nDim, input );
		fread( refinement_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&refinement_volume_min[i], sizeof(float) );
				reorder( (char *)&refinement_volume_max[i], sizeof(float) );
			}
		}

#ifdef STARFORM
		/* star formation volume */
		fread( &size, sizeof(int), 1, input );
		fread( star_formation_volume_min, sizeof(float), nDim, input );
		fread( star_formation_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&star_formation_volume_min[i], sizeof(float) );
				reorder( (char *)&star_formation_volume_max[i], sizeof(float) );
			}
		}
#endif /* STARFORM */
		
		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
		fread( &ncell0, sizeof(int), 1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&ncell0, sizeof(int) );
		}

		if ( ncell0 != (num_grid*num_grid*num_grid) ) {
			cart_error("File has different num_grid than compiled program");
		}
	}

	/* send header information to all other processors */
	MPI_Bcast( &endian, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &minlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &maxlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &step, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &box_size, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

#ifdef COSMOLOGY
	MPI_Bcast( &auni_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( (char *)&temp_cosmo, sizeof(struct CosmologyParameters), MPI_BYTE, MASTER_NODE, MPI_COMM_WORLD );
	cosmology_copy(&temp_cosmo);
#else
	MPI_Bcast( &h100, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &OmM0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

	/* h100 == 0 specifies user-defined units, > 0 implies ART standard */
	if(h100 > 0.0) {
		units_set_art(OmM0,h100,box_size);                                                                      
	} else {
		MPI_Bcast( &OmB0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &OmL0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

		units_set(OmM0,OmB0,OmL0);
	}
#endif /* COSMOLOGY */

	MPI_Bcast( tl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl_old, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
#ifdef HYDRO
	MPI_Bcast( level_sweep_dir, max_level-min_level+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
#endif /* HYDRO */

	MPI_Bcast( refinement_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( refinement_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );

#ifdef STARFORM
	MPI_Bcast( star_formation_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( star_formation_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
#endif /* STARFORM */

	if ( local_proc_id == file_parent ) {
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&size, sizeof(int) );
		}
		
		local_file_root_cells = size / sizeof(int);
		cellrefinedbuffer = cart_alloc(int, local_file_root_cells );

		num_read = fread( cellrefinedbuffer, sizeof(int), local_file_root_cells, input );

		if ( num_read != local_file_root_cells ) {
			cart_error("I/O error in read_grid_binary: num_read = %d, local_file_root_cells = %d",
				num_read, local_file_root_cells );
		}

		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < local_file_root_cells; i++ ) {
				reorder( (char *)&cellrefinedbuffer[i], sizeof(int) );
			}
		}

		if ( local_proc_id == MASTER_NODE ) {
			cell_counts = local_file_root_cells;

			file_sfc_index[0] = 0;
			for ( i = 1; i < num_output_files; i++ ) {
				/* need to know how many root cells are in each file */
				/* receive cellrefined */
				file_sfc_index[i] = cell_counts;
				MPI_Recv( &file_sfc_index[i+1], 1, MPI_INT, file_parent_proc[i],
					0, MPI_COMM_WORLD, &status );	
				cell_counts += file_sfc_index[i+1];
			}
			file_sfc_index[num_output_files] = cell_counts;

			cart_assert( cell_counts == num_root_cells );
		} else {
			/* send cellrefined array to MASTER_NODE */
			MPI_Send( &local_file_root_cells, 1, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD );
		}
	}

	/* send block information */
	MPI_Bcast( file_sfc_index, num_output_files+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	/* determine how many cells to expect from each processor */
	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_num_cells[proc] = 0;
	}

	for ( i = 0; i < num_output_files; i++ ) {
		if ( proc_sfc_index[local_proc_id] < file_sfc_index[i+1] &&
				proc_sfc_index[local_proc_id+1] >= file_sfc_index[i] ) {
			proc_num_cells[ file_parent_proc[i] ] = 
					min( proc_sfc_index[local_proc_id+1], file_sfc_index[i+1] ) - 
					max( proc_sfc_index[local_proc_id], file_sfc_index[i] );
		}
	}

	order = cart_alloc(int, num_cells_per_level[min_level] );
	cellrefined[local_proc_id] = cart_alloc(int, num_cells_per_level[min_level] );

	num_requests = 0;

	/* set up non-blocking receives */
	count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc_num_cells[proc] > 0 ) { 
			if ( proc == local_proc_id ) {
				/* copy data directly */
				start = max( 0, proc_sfc_index[local_proc_id] - file_sfc_index[file_index] );
				cart_assert( start >= 0 && start <= local_file_root_cells - proc_num_cells[local_proc_id] );
				for ( i = 0; i < proc_num_cells[local_proc_id]; i++ ) {
					cellrefined[local_proc_id][count+i] = cellrefinedbuffer[start+i];
				}
			} else {
				MPI_Irecv( &cellrefined[local_proc_id][count], proc_num_cells[proc], MPI_INT, proc,
					proc_num_cells[proc], MPI_COMM_WORLD, &requests[num_requests++] );
			}

			count += proc_num_cells[proc];
		}
	}

	cart_assert( count == num_cells_per_level[min_level] );

	/* need to avoid deadlocks! senders may also be receivers! */
	if ( local_proc_id == file_parent ) {
		start = 0;
		next_proc_index[proc] = 0;

		/* send all necessary information */
		for ( proc = 0; proc < num_procs; proc++ ) {
			next_proc_index[proc+1] = 0;

			if ( proc == local_proc_id ) { 
				for ( i = 0; i < proc_num_cells[local_proc_id]; i++ ) {
					if ( cellrefinedbuffer[start+i] > 1 ) {
						next_proc_index[proc+1]++;
					}
				}
				start += proc_num_cells[local_proc_id];
			} else if ( start < local_file_root_cells && proc != local_proc_id && 
					file_sfc_index[file_index]+start < proc_sfc_index[proc+1] &&
					file_sfc_index[file_index]+start >= proc_sfc_index[proc] ) {
				count = min( proc_sfc_index[proc+1], file_sfc_index[file_index+1] ) -
						max( file_sfc_index[file_index] + start, proc_sfc_index[proc] ); 
				cart_assert( count > 0 && start + count <= local_file_root_cells );
				for ( i = 0; i < count; i++ ) {
					if ( cellrefinedbuffer[start+i] > 1 ) {
						next_proc_index[proc+1]++;
					}
				}
				MPI_Isend( &cellrefinedbuffer[start], count, MPI_INT, proc, 
					count, MPI_COMM_WORLD, &requests[num_requests++] );
				start += count;
			}
		}

		cart_assert( start == local_file_root_cells );
	}

	/* wait for sends/receives to complete */
	MPI_Waitall( num_requests, requests, MPI_STATUSES_IGNORE );

	if ( local_proc_id == file_parent ) {
		cart_free( cellrefinedbuffer );
	}

	current_level_count = 0;
	next_level_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_cell_index[proc] = proc_sfc_index[local_proc_id] + current_level_count;
		proc_next_level_octs[proc] = 0;
	
		for ( i = 0; i < proc_num_cells[proc]; i++ ) {
			if ( cellrefined[local_proc_id][current_level_count] > 1 ) {
				ret = split_cell(current_level_count);
				if ( ret ) {
					cart_error("Unable to finish splitting root cells, ran out of octs?");
				}
				cart_assert( next_level_count < num_cells_per_level[min_level] );
				order[next_level_count++] = cell_child_oct[current_level_count];
				proc_next_level_octs[proc]++;
			}
			current_level_count++;
		}
	}

	cart_assert( next_level_count <= num_cells_per_level[min_level] );

	cart_free( cellrefined[local_proc_id] );
	//snl1
	read_grid_binary_top_level_vars(num_in_hydro_vars,skip_hvar,num_hydro_vars,hydro_vars,input,endian,file_parent,file_index,
			local_file_root_cells,page_size,proc_num_cells,proc_cell_index,file_sfc_index);
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
	read_grid_binary_top_level_vars(num_in_other_vars,skip_ovar,num_other_vars,other_vars,input,endian,file_parent,file_index,
			local_file_root_cells,page_size,proc_num_cells,proc_cell_index,file_sfc_index);
#endif /* GRAVITY || RADIATIVE_TRANSFER */

	/* now read levels */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		num_requests = 0;
		count = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_num_cells[proc] = proc_next_level_octs[proc]*num_children;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				MPI_Irecv( &first_page_count[proc], 1, MPI_LONG, proc,
						0, MPI_COMM_WORLD, &requests[num_requests++] );
			} else {
				first_page_count[proc] = 0;
			}

			proc_cell_index[proc] = count;
			count += proc_next_level_octs[proc];
			proc_next_level_octs[proc] = 0;
		}

		if ( local_proc_id == file_parent ) {
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}

			if ( size == sizeof(int) ) {
				fread( &size, sizeof(int), 1, input );

				if ( endian ) {
					reorder( (char *)&size, sizeof(int) );
				}

				total_cells[level] = (long)size;
			} else if ( size == sizeof(long) )  {
				fread( &total_cells[level], sizeof(long), 1, input );
				if ( endian ) {
					reorder( (char *)&total_cells[level], sizeof(long) );
				}
			} else {
				cart_error("File format error in %s, size = %d!", filename, size );
			}
			fread( &size, sizeof(int), 1, input );
			fread( &size, sizeof(int), 1, input );

			cart_debug("total_cells[%u] = %d", level, total_cells[level] );

			proc_first_index[0] = 0;
			for ( proc = 1; proc <= num_procs; proc++ ) {
				proc_first_index[proc] = (long)num_children*next_proc_index[proc] + 
								proc_first_index[proc-1];
				next_proc_index[proc] = 0;
			}

			cart_assert( proc_first_index[num_procs] == total_cells[level] );

			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( proc_first_index[proc+1]-proc_first_index[proc] > 0 && proc != local_proc_id ) {
					MPI_Isend( &proc_first_index[proc], 1, MPI_LONG, proc,
						0, MPI_COMM_WORLD, &requests[num_requests++] );
				}
			}

			cellrefinedbuffer = cart_alloc(int, min( total_cells[level], page_size ) );
		}

		MPI_Waitall( num_requests, requests, MPI_STATUS_IGNORE );

		num_requests = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_cur_cells[proc] = 0;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				/* set up receive */
				proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, 
								proc_num_cells[proc] );
				cellrefined[proc] = cart_alloc(int, min( page_size, proc_num_cells[proc] ) );
				MPI_Irecv( cellrefined[proc], proc_page_count[proc], MPI_INT, proc,
					proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				requests[proc] = MPI_REQUEST_NULL;
			}
		}

		current_level = order;
		order = cart_alloc(int, num_children*next_level_count );

		current_level_count = 0;
		current_read_count = 0;
		next_level_count = 0;
		continue_reading = 1;
		ready_to_read = 1;
	
		while ( continue_reading ) {
			flag = 0;

			if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read ) {
				/* read in a page */
				page_count = min( page_size, total_cells[level] - current_read_count );
				num_read = fread( cellrefinedbuffer, sizeof(int), page_count, input );

				if ( num_read != page_count ) {
					cart_error("I/O Error in read_grid_binary: num_read = %u, page_count = %u",
						num_read, page_count );
				}

				if ( endian ) {
					for ( i = 0; i < page_count; i++ ) {
						reorder( (char *)&cellrefinedbuffer[i], sizeof(int) );
					}
				}

				/* send info to other processors */
				start = 0;
				num_send_requests = 0;

				for ( proc = 0; proc < num_procs; proc++ ) {
					if ( start < page_count && current_read_count + start >= proc_first_index[proc] &&
							current_read_count + start < proc_first_index[proc+1] ) {
						
						count = min( proc_first_index[proc+1], current_read_count+page_count ) -
							max( proc_first_index[proc], current_read_count + start );

						cart_assert( count > 0 && count <= page_count );

						for ( i = 0; i < count; i++ ) {
							cart_assert( start + i < page_count );
							if ( cellrefinedbuffer[start+i] ) {
								next_proc_index[proc+1]++;
							}
						}

						if ( proc == local_proc_id ) {
							/* copy into local buffer */
							proc_page_count[local_proc_id] = count;
							cellrefined[local_proc_id] =
								cart_alloc(int, proc_page_count[local_proc_id] );
							for ( i = 0; i < proc_page_count[local_proc_id]; i++ ) {
								cart_assert( start + i < page_count );
								cellrefined[local_proc_id][i] = cellrefinedbuffer[start+i];
							}
						} else {
							MPI_Isend( &cellrefinedbuffer[start], count, MPI_INT, proc,
								count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
						}

						start += count;
					}
				}

				cart_assert( start == page_count );

				current_read_count += page_count;
				ready_to_read = 0;

				/* split local cells first */
				if ( proc_page_count[local_proc_id] > 0 ) {
					flag = 1;
					proc = local_proc_id;
				} else if ( num_requests > 0 ) {
					/* see if we've received anything */
					MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
				}

			} else if ( num_requests > 0 ) {
				/* wait for a receive to complete */
				MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
				num_requests--;
				flag = 1;
			}

			if ( flag == 1 && proc != MPI_UNDEFINED ) {
				cart_assert( proc >= 0 && proc < num_procs );

				/* split cells in received page */
				i = 0;
				while ( i < proc_page_count[proc] ) {
					cart_assert( proc_cell_index[proc] + proc_cur_cells[proc]/num_children 
							< num_cells_per_level[level]/num_children );
					ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children ]; 

					cart_assert( ioct >= 0 && ioct < num_octs );
					cart_assert( oct_level[ioct] == level );

					for ( j = 0; j < num_children; j++ ) {
						icell = oct_child( ioct, j );
						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						if ( cellrefined[proc][i] == 1 ) {
							ret = split_cell(icell);

							if ( ret ) {
								cart_error("Unable to finish splitting cells, ran out of octs?");
							}

							order[ proc_cell_index[proc]*num_children +
								proc_next_level_octs[proc] ] = cell_child_oct[icell];
							next_level_count++;
							proc_next_level_octs[proc]++;

						}

						i++;
						proc_cur_cells[proc]++;
						current_level_count++;
					}
				}

				/* if necessary, start a new receive from that proc */
				if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
					proc_page_count[proc] = min( page_size, 
						proc_num_cells[proc] - proc_cur_cells[proc] );
					cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
					MPI_Irecv( cellrefined[proc], proc_page_count[proc], MPI_INT,
						proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
					num_requests++;
				} else {
					proc_page_count[proc] = 0;
					cart_free( cellrefined[proc] );
					requests[proc] = MPI_REQUEST_NULL;
				}

				flag = 0;
			}

			/* check if we're done reading */
			if ( local_proc_id == file_parent ) {
				if ( num_requests > 0 ) {
					/* see if sends have completed */
					MPI_Testall( num_send_requests, send_requests, 
							&ready_to_read, MPI_STATUSES_IGNORE );
				} else {
					MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
					num_send_requests = 0;
					ready_to_read = 1;
				}

				if ( current_read_count >= total_cells[level] &&
						current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			} else {
				if ( current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( !ready_to_read ) {
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
			}

			cart_free( cellrefinedbuffer );

			fread( &size, sizeof(int), 1, input );
		}

		/* repack the order array */
		count = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			for ( i = 0; i < proc_next_level_octs[proc]; i++ ) {
				cart_assert( proc_cell_index[proc]*num_children + i  >= count );
				cart_assert( order[ proc_cell_index[proc]*num_children + i ] >= 0 &&
						order[ proc_cell_index[proc]*num_children + i ] < num_octs );
				cart_assert( oct_level[ order[ proc_cell_index[proc]*num_children + i ] ] == level+1 );

				order[count++] = order[ proc_cell_index[proc]*num_children + i ];
			}
		}

		cart_assert( count == next_level_count );

		read_grid_binary_lower_level_vars(num_in_hydro_vars,skip_hvar,num_hydro_vars,hydro_vars,input,endian,file_parent,
					file_index,total_cells,page_size,proc_num_cells,level,first_page_count,
					proc_first_index,proc_cell_index,current_level);

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		read_grid_binary_lower_level_vars(num_in_other_vars,skip_ovar,num_other_vars,other_vars,input,endian,file_parent,file_index,
					total_cells,page_size,proc_num_cells,level,first_page_count,proc_first_index,
					proc_cell_index,current_level);
#endif /* GRAVITY || RADIATIVE_TRANSFER */

		cart_free( current_level );		
	}

	if ( local_proc_id == file_parent ) {
		fclose( input );
	}

	if(other_vars != 0) cart_free(other_vars);

	cart_free( order );

#ifdef RADIATIVE_TRANSFER
	/* Load RF data */
	rtReadRadiationFieldData(filename,1);
#endif
}

#ifdef BLASTWAVE_FEEDBACK
extern double blastwave_time_floor;
extern double blastwave_time_cut;
extern double blastwave_time_code_max;
#endif /* BLASTWAVE_FEEDBACK */

float set_skipped_var(int jskip){
#ifdef BLASTWAVE_FEEDBACK
  if(jskip == HVAR_BLASTWAVE_TIME - HVAR_GAS_DENSITY){ 
    return (float)1.1*blastwave_time_floor;
  }else{
    cart_error("skipped variables %d not set",jskip);
  }
#else
  cart_error("skipped variables %d not set",jskip);
  return -1;
#endif /* BLASTWAVE_FEEDBACK */
}

void read_grid_binary_top_level_vars(int num_in_vars, int jskip, int num_out_vars, int *out_var, FILE *input, int endian, 
		int file_parent, int file_index, int local_file_root_cells, int page_size, int *proc_num_cells, 
		long *proc_cell_index, int *file_sfc_index) {
        int i, j, m, jout;
	int size;
	float *cellvars[MAX_PROCS], *cellvars_buffer, *cellvars_read_buffer;
	int num_requests, num_send_requests;
	int proc;
	int proc_cur_cells[MAX_PROCS];
	int proc_page_count[MAX_PROCS];
	long current_level_count, current_read_count;
	int continue_reading, ready_to_read;
	int num_read, flag;
	int page_count;
	long start, count;
	int icell;
	MPI_Request requests[2*MAX_PROCS];
	MPI_Request send_requests[MAX_PROCS];

	if(num_out_vars < 1) return;
	cart_assert( num_in_vars >= num_out_vars || (jskip >= 0 && (num_in_vars >= num_out_vars-1)) );

	if ( local_proc_id == file_parent ) {
		fread( &size, sizeof(int), 1, input );
		cellvars_buffer = cart_alloc(float, num_out_vars*min( local_file_root_cells, page_size ) );
		cellvars_read_buffer = cart_alloc(float, num_in_vars*min( local_file_root_cells, page_size ) );
	}

	num_requests = 0;

	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_cur_cells[proc] = 0;

		if ( proc_num_cells[proc] > 0 && proc != local_proc_id )  {
			/* set up receive */
			proc_page_count[proc] = min( proc_num_cells[proc], 
					page_size - ( proc_cell_index[proc] - 
					file_sfc_index[(int)((proc * num_output_files)/ num_procs)]) % page_size);
			cellvars[proc] = cart_alloc(float, num_out_vars*min( page_size, proc_num_cells[proc] ) );
			MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], 
					MPI_COMM_WORLD, &requests[proc] );
			num_requests++;
		} else {
			proc_page_count[proc] = 0;
			requests[proc] = MPI_REQUEST_NULL;
		}
	}

	current_level_count = 0;
	current_read_count = 0;
	continue_reading = 1;
	ready_to_read = 1;

	while ( continue_reading ) {
		flag = 0;

		if ( local_proc_id == file_parent && current_read_count < local_file_root_cells && ready_to_read ) {
			/* read in a page */
			page_count = min( page_size, local_file_root_cells - current_read_count );
			num_read = fread( cellvars_read_buffer, sizeof(float), num_in_vars*page_count, input );

			if ( num_read != num_in_vars*page_count ) {
				cart_error("I/O Error in read_grid_binary: num_read = %u, num_in_vars*page_count = %u", 
						num_read, num_in_vars*page_count );
			}

			for(i=0; i<page_count; i++)
			  {
			    jout = 0; for(j=0; j<max(num_in_vars,num_out_vars); j++)
			      {
				while( jout == jskip ){
				  cellvars_buffer[num_out_vars*i+jskip] = set_skipped_var(jskip);
				  jout++;
				}
				if( jout < num_out_vars ) {cellvars_buffer[num_out_vars*i+jout] = cellvars_read_buffer[num_in_vars*i+j];}
				jout++;
			      }
			  }

			if ( endian ) {
				for ( i = 0; i < num_out_vars*page_count; i++ ) {
					reorder( (char *)&cellvars_buffer[i], sizeof(float) );
				}
			}

			/* send info to other processors */
			start = 0;
			num_send_requests = 0;
			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( start < page_count && 
						file_sfc_index[file_index]+current_read_count+start < proc_sfc_index[proc+1] &&
						file_sfc_index[file_index]+current_read_count+start >= proc_sfc_index[proc] ) {
		  
					count = min( proc_sfc_index[proc+1], file_sfc_index[file_index]+current_read_count+page_count ) - 
							max( proc_sfc_index[proc],  file_sfc_index[file_index] + current_read_count + start );
					cart_assert( count > 0 && count <= page_count );

					if ( proc == local_proc_id ) {
						/* copy into local buffer */
						proc_page_count[local_proc_id] = count;
						cellvars[local_proc_id] = cart_alloc(float, num_out_vars*proc_page_count[local_proc_id] );

						for ( i = 0; i < num_out_vars*proc_page_count[local_proc_id]; i++ ) {
							cart_assert( num_out_vars*start + i < num_out_vars*page_count );
							cellvars[local_proc_id][i] = cellvars_buffer[num_out_vars*start+i];
						}
					} else { 
						MPI_Isend( &cellvars_buffer[num_out_vars*start], num_out_vars*count, MPI_FLOAT, proc, count, 
								MPI_COMM_WORLD, &send_requests[num_send_requests++] );
					}
		  
					start += count;
				}
			}

			cart_assert( start == page_count );
			current_read_count += page_count;
			ready_to_read = 0;

			/* unpack local cells first */
			if ( proc_page_count[local_proc_id] > 0 ) {
				flag = 1;
				proc = local_proc_id;
			} else if ( num_requests > 0 ) {
				/* see if we've received anything */
				MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
			}
		} else if ( num_requests > 0 ) {
			/* wait for a receive to complete */
			MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
			num_requests--;
			flag = 1;
		}

		if ( flag == 1 && proc != MPI_UNDEFINED ) {
			cart_assert( proc >= 0 && proc < num_procs );
			
			/* unpack received page */
			i = 0;
			while ( i < num_out_vars*proc_page_count[proc] ) {
				icell = root_cell_location( proc_cell_index[proc] + proc_cur_cells[proc]);
				cart_assert( icell >= 0 && icell < num_cells_per_level[min_level] );
	      
				for ( m = 0; m < num_out_vars; m++ ) {
					cell_var(icell,out_var[m]) = cellvars[proc][i++];
				}

				current_level_count++;
				proc_cur_cells[proc]++;
			}

			/* if necessary, start a new receive for that proc */
			if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
				proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
				cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
				MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], 
						MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				cart_free( cellvars[proc] );
				requests[proc] = MPI_REQUEST_NULL;
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( num_requests > 0 ) {
				/* see if sends have completed, if they have we can continue to read */
				MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
			} else {
				/* not waiting on any receives, so we can't do anything until
				 * we read in additional data */
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
				num_send_requests = 0;
				ready_to_read = 1;
			}

			if ( current_read_count >= local_file_root_cells && current_level_count >= num_cells_per_level[min_level] ) {
				continue_reading = 0;
			}
		} else {
			if ( current_level_count >= num_cells_per_level[min_level] ) {
				continue_reading = 0;
			}
		}
	}

	if ( local_proc_id == file_parent ) {
		if ( !ready_to_read ) {
			MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
		}

		cart_free( cellvars_read_buffer );
		cart_free( cellvars_buffer );
		fread( &size, sizeof(int), 1, input );
	}
}


void read_grid_binary_lower_level_vars(int num_in_vars, int jskip, int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, long *total_cells, int page_size, int *proc_num_cells, int level, long *first_page_count, long *proc_first_index, long *proc_cell_index, int *current_level)
{
  int i, j, m, jout;
  int size;
  float *cellvars[MAX_PROCS], *cellvars_buffer, *cellvars_read_buffer;
  int num_requests, num_send_requests;
  int proc;
  int proc_cur_cells[MAX_PROCS];
  int proc_page_count[MAX_PROCS];
  long current_level_count, current_read_count;
  int continue_reading, ready_to_read;
  int num_read, flag;
  int page_count;
  long start, count;
  int icell, ioct;
  MPI_Request requests[2*MAX_PROCS];
  MPI_Request send_requests[MAX_PROCS];


  if(num_out_vars < 1) return;
  cart_assert( num_in_vars >= num_out_vars || (jskip >= 0 && (num_in_vars >= num_out_vars-1)) );

  if ( local_proc_id == file_parent )
    {
      fread( &size, sizeof(int), 1, input );
      cellvars_buffer = cart_alloc(float, num_out_vars*min( total_cells[level], page_size ) );
      cellvars_read_buffer = cart_alloc(float, num_in_vars*min( total_cells[level], page_size ) );
    }
	
  num_requests = 0;

  for ( proc = 0; proc < num_procs; proc++ )
    {
      proc_cur_cells[proc] = 0;

      if ( proc_num_cells[proc] > 0 && proc != local_proc_id )
	{
	  /* set up receive */
	  proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, proc_num_cells[proc] );
	  cellvars[proc] = cart_alloc(float, num_out_vars*min( page_size, proc_num_cells[proc] ) );
	  MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
	  num_requests++;
	}
      else
	{
	  proc_page_count[proc] = 0;
	  requests[proc] = MPI_REQUEST_NULL;
	}
    }
	
  current_level_count = 0;
  current_read_count = 0;
  continue_reading = 1;
  ready_to_read = 1;

  while ( continue_reading )
    {
      flag = 0;

      if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read )
	{
	  /* read in a page */
	  page_count = min( page_size, total_cells[level] - current_read_count );
	  num_read = fread( cellvars_read_buffer, sizeof(float), num_in_vars*page_count, input );
	
	  if ( num_read != num_in_vars*page_count )
	    {
	      cart_error("I/O Error in read_grid_binary: num_read = %u, num_out_vars*page_count = %u", num_read, num_out_vars*page_count );
	    }
	
	  for(i=0; i<page_count; i++)
	    {
	      jout=0; for(j=0; j<max(num_in_vars,num_out_vars); j++)
		{
		  while( jout == jskip ){
		    cellvars_buffer[num_out_vars*i+jskip] = set_skipped_var(jskip);
		    jout++;
		  }
		  if( jout < num_out_vars ) {cellvars_buffer[num_out_vars*i+jout] = cellvars_read_buffer[num_in_vars*i+j];}
		  jout++;
		}
	    }

	  if ( endian )
	    {
	      for ( i = 0; i < num_out_vars*page_count; i++ )
		{
		  reorder( (char *)&cellvars_buffer[i], sizeof(float) );
		}
	    }
	
	  /* send info to other processors */
	  start = 0;
	  num_send_requests = 0;
	  for ( proc = 0; proc < num_procs; proc++ )
	    {
	      if ( start < page_count && current_read_count + start >= proc_first_index[proc] && current_read_count + start < proc_first_index[proc+1] )
		{
				
		  count = min( proc_first_index[proc+1], current_read_count+page_count ) - max( proc_first_index[proc], current_read_count + start );
		  cart_assert( count > 0 && count <= page_count );

		  if ( proc == local_proc_id )
		    {
		      /* copy into local buffer */
		      proc_page_count[local_proc_id] = count;

		      cellvars[local_proc_id] = cart_alloc(float, num_out_vars*proc_page_count[local_proc_id] );
	
		      for ( i = 0; i < num_out_vars*proc_page_count[local_proc_id]; i++ )
			{
			  cart_assert( num_out_vars*start + i < num_out_vars*page_count );
			  cellvars[local_proc_id][i] = cellvars_buffer[num_out_vars*start+i];
			}
		    }
		  else
		    {
		      MPI_Isend( &cellvars_buffer[num_out_vars*start], num_out_vars*count, MPI_FLOAT, proc, count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
		    }
	
		  start += count;
		}
	    }
	
	  cart_assert( start == page_count );
	  current_read_count += page_count;
	  ready_to_read = 0;
	
	  /* unpack local cells first */
	  if ( proc_page_count[local_proc_id] > 0 )
	    {
	      flag = 1;
	      proc = local_proc_id;
	    }
	  else if ( num_requests > 0 )
	    {
	      /* see if we've received anything */
	      MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
	    }
	}
      else if ( num_requests > 0 )
	{
	  /* wait for a receive to complete */
	  MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
	  num_requests--;
	  flag = 1;
	}
	
      if ( flag == 1 && proc != MPI_UNDEFINED )
	{
	  cart_assert( proc >= 0 && proc < num_procs );

	  /* unpack received page */
	  i = 0;
	
	  while ( i < num_out_vars*proc_page_count[proc] )
	    {
	      ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children];
	      cart_assert( ioct >= 0 && ioct < num_octs );

	      for ( j = 0; j < num_children; j++ )
		{	
		  icell = oct_child( ioct, j );
		  cart_assert( icell >= 0 && icell < num_cells );
		  cart_assert( cell_level(icell) == level );

		  for ( m = 0; m < num_out_vars; m++ )
		    {
		      cell_var(icell,out_var[m]) = cellvars[proc][i++];
		    }

		  
		  current_level_count++;
		  proc_cur_cells[proc]++;
		}
	    }

	  /* if necessary, start a new receive for that proc */
	  if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id )
	    {
	      proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
	      cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
	      MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
	      num_requests++;
	    }
	  else
	    {
	      proc_page_count[proc] = 0;
	      cart_free( cellvars[proc] );
	      requests[proc] = MPI_REQUEST_NULL;
	    }
	
	  flag = 0;
	}

      if ( local_proc_id == file_parent )
	{
	  if ( num_requests > 0 )
	    {
	      /* see if sends have completed */
	      MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
	    }
	  else
	    {
	      MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	      num_send_requests = 0;
	      ready_to_read = 1;
	    }
	  
	  if ( current_read_count >= total_cells[level] && current_level_count >= num_cells_per_level[level] )
	    {
	      continue_reading = 0;
	    }
	}
      else
	{
	  if ( current_level_count >= num_cells_per_level[level] )
	    {
	      continue_reading = 0;
	    }
	}
    }

  if ( local_proc_id == file_parent )
    {
      if ( !ready_to_read )
	{
	  MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	}
	
      cart_free( cellvars_buffer );
      fread( &size, sizeof(int), 1, input );
    }
}

