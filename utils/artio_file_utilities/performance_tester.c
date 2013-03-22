#include "config.h"


#include <mpi.h>
#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "control_parameter.h"
#include "io.h"
#include "io_artio.h"


extern const char* executable_name;
extern char output_directory_d[]; 
extern char logfile_directory_d[]; 


char apt_output_directory[CONTROL_PARAMETER_STRING_LENGTH] = "";
const char *apt_label = NULL;


void init()
{
  const char *str;
  char c, *tmp;
  int n;
  float v;

  if((!is_option_present("job-name","j",1) && !is_option_present("input","i",1)))
    {
      cart_error("Usage: %s -i/--input=<filename> [options].\n"
		 "   or: %s -j/--job-name=<name> [options].\n"
		 "Valid options:\n"
		 "  -j, --job-name=<name>              set the job name for the fileset\n" 
		 "  -d, --data-directory=<dir>         set the input data directory where files are\n" 
		 "                                     located (default is the current directory)\n"
		 "  -o, --output-directory=<dir>       set the data directory where converted output\n"
		 "                                     files are written (default is data-directory).\n"
		 "  -l, --file-label=<label>           set the label for the fileset\n"
		 "  -i, --input=<filename>             set the input grid filename; this option\n"
		 "                                     combines -j, -d, and -l in one option, i.e.\n"
		 "                                     <filename> = <dir>/<name>_<label>.art.\n"
                 "  -na, --num-artio_files=<number>    set the number of artio files (output only).\n"
		 "  ...                                other options (TBD)\n"
		,executable_name,executable_name);
    }

  control_parameter_add2(control_parameter_string,(char *)jobname,"job-name","jobname","");

  str = extract_option1("input","i",NULL);
  if(str != NULL)
    {
      tmp = strrchr(str,'.');
      if(tmp != NULL)
	{
	  /*
	  //  Chop off the extension
	  */
	  *tmp = 0;
	  tmp++;
	}
      if(tmp==NULL || strcmp(tmp,"art")!=0)
	{
	  cart_error("Input ARTIO file must end with .art.");
	}

      tmp = strrchr(str,'/');
      if(tmp != NULL)
	{
	  /*
	  //  Detach the leading directory path
	  */
	  *tmp = 0;
	  strncpy(output_directory_d,str,CONTROL_PARAMETER_STRING_LENGTH);
	  output_directory_d[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
	  str = tmp + 1;
	  if(strlen(str) == 0)
	    {
	      cart_error("Missing the file name after the last /.");
	    }
	}

      tmp = strrchr(str,'_');
      if(tmp != NULL)
	{
	  if(sscanf(tmp,"_a%f%c",&v,&c)==1 || sscanf(tmp,"_%d%c",&n,&c)==1)
	    {
	      apt_label = tmp + 1;
	      *tmp = 0;
	    }
	}

      set_jobname(str);
    }

  str = extract_option1("job-name","j",NULL);
  if(str != NULL)
    {
      set_jobname(str);
    }

  str = extract_option1("data-directory","d",NULL);
  if(str != NULL)
    {
      strncpy(output_directory_d,str,CONTROL_PARAMETER_STRING_LENGTH-1);
      output_directory_d[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
    }

  str = extract_option1("output-directory","o",NULL);
  if(str != NULL)
    {
      strncpy(apt_output_directory,str,CONTROL_PARAMETER_STRING_LENGTH-1);
      apt_output_directory[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
    }

  str = extract_option1("file-label","l",NULL);
  if(str != NULL)
    {
      apt_label = str;
    }
}


void read_file()
{
  double wtime = MPI_Wtime();

  read_artio_restart(apt_label);

  cart_debug("Reading time: %lg sec",MPI_Wtime()-wtime);
}


void write_file()
{
  const char *str;
  double wtime;
  char c;
  int n;

  if ( strcmp(apt_output_directory,"") != 0 ) {
      cart_debug("Setting output directory to %s", apt_output_directory);
      strncpy(output_directory_d,apt_output_directory,CONTROL_PARAMETER_STRING_LENGTH-1);
      output_directory_d[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
  }

  while((str = extract_option1("num-artio-files","na",NULL)) != NULL)
    {
      if(sscanf(str,"%d%c",&n,&c)!=1 || n<1)
	{
	  cart_error("Option --num-artio-files must have a positive integer number as its argument.");
	}
      num_artio_grid_files = n;
#ifdef PARTICLES
      num_artio_particle_files = n;
#endif /* PARTICLES */

      wtime = MPI_Wtime();
      write_artio_restart(WRITE_GENERIC,WRITE_GENERIC,WRITE_GENERIC);
      cart_debug("Writing time: %lg sec",MPI_Wtime()-wtime);
    }

}
