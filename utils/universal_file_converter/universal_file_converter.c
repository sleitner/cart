#include "config.h"


#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "control_parameter.h"
#include "io.h"
#include "io_artio.h"
#include "io_cart.h"
#include "parallel.h"
#include "../extra/hart_io.h"


extern const char* executable_name;
extern char output_directory_d[]; 
extern char logfile_directory_d[]; 

char ufc_output_directory[CONTROL_PARAMETER_STRING_LENGTH] = "";

#define H2C   1
#define H2A   2
#define C2A   3


int ufc_mode = 0;
const char *ufc_label = NULL;
const char *ufc_hart_fname = NULL;
int ufc_flag = 3;     /* convert everything */


void init()
{
  const char *str, *outdir;
  char c, *tmp;
  int n, mode;
  float v;

  if(strcmp(options[0],"hart-to-cart")==0 || strcmp(options[0],"h2c")==0)
    {
      ufc_mode = +H2C;
    }

  if(strcmp(options[0],"cart-to-hart")==0 || strcmp(options[0],"c2h")==0)
    {
      if(ufc_mode != 0) cart_error("Only one conversion specification is allowed.");
      ufc_mode = -H2C;
    }

  if(strcmp(options[0],"hart-to-artio")==0 || strcmp(options[0],"h2a")==0)
    {
      if(ufc_mode != 0) cart_error("Only one conversion specification is allowed.");
      ufc_mode = +H2A;
    }

  if(strcmp(options[0],"artio-to-hart")==0 || strcmp(options[0],"a2h")==0)
    {
      if(ufc_mode != 0) cart_error("Only one conversion specification is allowed.");
      ufc_mode = -H2A;
    }

  if(strcmp(options[0],"cart-to-artio")==0 || strcmp(options[0],"c2a")==0)
    {
      if(ufc_mode != 0) cart_error("Only one conversion specification is allowed.");
      ufc_mode = +C2A;
    }

  if(strcmp(options[0],"artio-to-cart")==0 || strcmp(options[0],"a2c")==0)
    {
      if(ufc_mode != 0) cart_error("Only one conversion specification is allowed.");
      ufc_mode = -C2A;
    }

  options++;
  num_options--;

  if(ufc_mode==0 || (!is_option_present("job-name","j",1) && !is_option_present("input","i",1)))
    {
      cart_error("Usage: %s <conversion-spec> -i/--input=<filename> [options].\n"
		 "   or: %s <conversion-spec> -j/--job-name=<name> [options].\n"
		 "Valid conversion specifications:\n"
		 "  hart-to-cart or cart-to-hart (shorthand h2c/c2h)\n"
		 "  hart-to-artio or artio-to-hart (shorthand h2a/a2h)\n"
		 "  cart-to-artio or artio-to-cart (shorthand c2a/a2c)\n"
		 "Valid options:\n"
		 "  -j, --job-name=<name>              set the job name for the fileset\n" 
		 "  -d, --data-directory=<dir>         set the input data directory where files are\n" 
		 "                                     located (default is the current directory)\n"
         "  -o, --output-directory=<dir>       set the data directory where converted output\n"
         "                                     files are written (default is data-directory).\n"
		 "  -l, --file-label=<label>           set the label for the fileset\n"
		 "  -i, --input=<filename>             set the input grid filename; this option\n"
		 "                                     combines -j, -d, and -l in one option, i.e.\n"
		 "                                     <filename> = <dir>/<name>_<label>.ext where\n"
		 "                                     .ext is either .d or .art\n" 
		 "  -hf, --hart-file-name=<name>       set the name for the HART grid file\n"
		 "                                     (default uses a .dh suffix)\n"
		 "  -nc, --num-cart_files=<number>     set the number of cart files (output only)\n"
		 "  -na, --num-artio_files=<number>    set the number of artio files (output only)\n"
#ifdef PARTICLES
		 "  -p, --particle-only                convert particle files only\n"
		 "                                     (for artio-to-* conversion)\n"
		 "  -g, --grid-only                    convert grid files only\n"
		 "                                     (for artio-to-* conversion)\n"
		 "  -nrow, --cart-num-row=<number>     set the NROW parameter of hart/cart\n"
		 "                                     particle files (output only)\n"
		 "  -pfm, --particle-file-mode=<mode>  set the old-style particle file mode\n"
		 "                                     for reading legacy files\n"
         "                                     0 = double-prec positions, double-prec times (default),\n"
         "                                     1 = double-prec positions, single-prec times,\n"
         "                                     2 = single-prec positions, single-prec times\n"
#endif
		 "  -gfm, --grid-file-mode=<mode>      set the old-style grid file mode for\n"
		 "                                     reading legacy files\n"
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
      if(tmp==NULL || (strcmp(tmp,"d")!=0 && strcmp(tmp,"art")!=0))
	{
	  cart_error("Input grid file must end on .d or .art.");
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
	      ufc_label = tmp + 1;
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
      strncpy(ufc_output_directory,str,CONTROL_PARAMETER_STRING_LENGTH-1);
      ufc_output_directory[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
    }

  str = extract_option1("file-label","l",NULL);
  if(str != NULL)
    {
      ufc_label = str;
    }

  ufc_hart_fname = extract_option1("hart-file-name","hf",NULL);
  if(ufc_hart_fname != NULL)
    {
      if(ufc_mode==C2A || ufc_mode==-C2A)
	{
	  cart_error("Option --hart-file-name should be used with HART file conversion only.");
	}
    }

  str = extract_option1("num-cart-files","nc",NULL);
  if(str != NULL)
    {
      if(sscanf(str,"%d%c",&n,&c)!=1 || n<1 || n>num_procs)
	{
	  cart_error("Option --num-cart-files must have an integer number between 1 and %d as its argument.",num_procs);
	}

      if ( ufc_mode == C2A || ufc_mode == -H2C ) {
		num_cart_input_files = n;
	  } else if ( ufc_mode == H2C || ufc_mode == -C2A ) {
        num_cart_output_files = n;
	  } else {
		cart_error("Option --num-cart-files should be used with CART file conversion only.");
      }
    }
  
  str = extract_option1("num-artio-files","na",NULL);
  if(str != NULL)
    {
      if(sscanf(str,"%d%c",&n,&c)!=1 || n<1)
	{
	  cart_error("Option --num-artio-files must have a positive integer number as its argument.");
	}
      num_artio_grid_files = n;
#ifdef PARTICLES
      num_artio_particle_files = n;
#endif /* PARTICLES */
    }

#ifdef PARTICLES
  str = extract_option1("cart-num-row","nrow",NULL);
  if(str != NULL)
    {
      if(sscanf(str,"%d%c",&n,&c)!=1 || n<1)
	{
	  cart_error("Option --cart-num-row must have a positive integer number as its argument.",num_procs);
	}
      cart_particle_num_row = n;      
    }
  
  str = extract_option0("particle-only","p");
  if(str != NULL)
    {
      if(ufc_mode!=-H2A && ufc_mode!=-C2A)
	{
	  cart_error("Option --particle-only should be used with artio-to-* conversion only.");
	}
      ufc_flag = 2;
    }

  str = extract_option0("grid-only","p");
  if(str != NULL)
    {
      if(ufc_mode!=-H2A && ufc_mode!=-C2A)
	{
	  cart_error("Option --grid-only should be used with artio-to-* conversion only.");
	}
      if(ufc_flag != 3)
	{
	  cart_error("Options --particle-only and --grid-only are mutually exclusive.");
	}
      ufc_flag = 1;
    }

  /*
  //  Support -pfm/--particle-file-mode=<mode> option
  */
  str = extract_option1("particle-file-mode","pfm",NULL);
  if(str != NULL)
    {
      if(sscanf(str,"%d",&mode)==0 || mode<0 || mode>2)
	{
	  cart_error("--particle-file-mode=<mode> option requires an integer <mode> between 0 and 2 as an argument.");
	}
      set_cart_particle_file_mode(mode);
    }
#endif /* PARTICLES */

  /*
  //  Support -gfm/--grid-file-mode=<mode> option
  */
  str = extract_option1("grid-file-mode","gfm",NULL);
  if(str != NULL)
    {
      if(sscanf(str,"%d",&mode)==0 || mode<0 || mode>4)
	{
	  cart_error("--grid-file-mode=<mode> requires an integer <mode> between 0 and 3 as an argument.");
	}
      set_cart_grid_file_mode(mode);
    }

  die_on_unknown_options();

  if(ufc_hart_fname == NULL)
    {
		if ( strcmp(ufc_output_directory,"") != 0 ) {
			outdir = ufc_output_directory;
		} else {
			outdir = output_directory;
		}
	
      n = strlen(outdir) + 1 + strlen(jobname);
      if(ufc_label != NULL) n += 1 + strlen(ufc_label) + 4;
      str = cart_alloc(char,n);

      if(ufc_label == NULL)
	{
	  sprintf((char *)str,"%s/%s.dh",outdir,jobname);
	}
      else
	{
	  sprintf((char *)str,"%s/%s_%s.dh",outdir,jobname,ufc_label);
	}
      ufc_hart_fname = str;
    }
}


void read_file()
{
  double wtime = MPI_Wtime();

  switch(ufc_mode)
    {
    case +H2C:
    case +H2A:
      {
	read_hart_grid_binary((char *)ufc_hart_fname);
	break;
      }
    case -H2C:
    case +C2A:
      {
	read_cart_restart(ufc_label);
	break;
      }
    case -H2A:
    case -C2A:
      {
	read_artio_restart(ufc_label);
	break;
      }
    }

  cart_debug("Reading time: %lg sec",MPI_Wtime()-wtime);
}


void write_file()
{
  double wtime = MPI_Wtime();

  if ( strcmp(ufc_output_directory,"") != 0 ) {
      cart_debug("Setting output directory to %s", ufc_output_directory);
      strncpy(output_directory_d,ufc_output_directory,CONTROL_PARAMETER_STRING_LENGTH-1);
      output_directory_d[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
  }

  switch(ufc_mode)
    {
    case -H2A:
      {
	if(ufc_flag & 2) write_cart_restart(NO_WRITE,WRITE_SAVE,NO_WRITE);
      }
    case -H2C:
      {
	write_hart_grid_binary((char *)ufc_hart_fname);
	break;
      }
    case +H2C:
    case -C2A:
      {
	write_cart_restart((ufc_flag & 1) ? WRITE_SAVE : NO_WRITE,(ufc_flag & 2) ? WRITE_SAVE : NO_WRITE,WRITE_SAVE);
	break;
      }
    case +H2A:
    case +C2A:
      {
	write_artio_restart(WRITE_SAVE,WRITE_SAVE,WRITE_SAVE);
	break;
      }
    }

  cart_debug("Writing time: %lg sec",MPI_Wtime()-wtime);
}
