#include "config.h"


#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "control_parameter.h"
#include "io.h"
#include "io_artio.h"
#include "io_cart.h"
#include "../extra/hart_io.h"


extern const char* executable_name;
extern char output_directory_d[]; 
extern char logfile_directory_d[]; 


#define H2C   1
#define H2A   2
#define C2A   3


int ufc_mode = 0;
const char *ufc_label = NULL;
const char *ufc_hart_fname = NULL;
int ufc_flag = 3;     /* convert everything */


void init()
{
  const char *str;
  char c;
  int n;

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

  if(ufc_mode==0 || !is_option_present("job-name","j",1))
    {
      cart_error("Usage: %s <conversion-spec> -j/--job-name=<name> [options].\n"
		 "Valid conversion specifications:\n"
		 "  hart-to-cart or cart-to-hart (shorthand h2c/c2h)\n"
		 "  hart-to-artio or artio-to-hart (shorthand h2a/a2h)\n"
		 "  cart-to-artio or artio-to-cart (shorthand c2a/a2c)\n"
		 "Valid options:\n"
		 "  -d, --data-directory=<dir>       set data directory where files are located (default is the current directory)\n"
		 "  -l, --file-label=<label>         set the label for the fileset\n"
		 "  -hf, --hart-file-name=<name>     set the name for the HART grid file (default uses a .dh suffix)"
		 "  -nc, --num-cart_files=<number>   set the number of cart files\n"
		 "  -na, --num-artio_files=<number>  set the number of artio files\n"
#ifdef PARTICLES
		 "  -p, --particle-only              convert particle files only (for artio-to-* conversion)\n"
		 "  -g, --grid-only                  convert grid files only (for artio-to-* conversion)\n"
		 "  -nrow, --cart-num-row=<number>   set the NROW parameter of hart/cart particle files\n"
#endif
,executable_name);
    }

  str = extract_option1("job-name","j",NULL);

  control_parameter_add2(control_parameter_string,(char *)jobname,"job-name","jobname","");
  set_jobname(str);

  str = extract_option1("data-directory","d",NULL);
  if(str != NULL)
    {
      strncpy(output_directory_d,str,CONTROL_PARAMETER_STRING_LENGTH);
      output_directory_d[CONTROL_PARAMETER_STRING_LENGTH-1] = 0;
    }

  ufc_label = extract_option1("file-label","l",NULL);

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
      num_cart_output_files = n;      
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
#endif

  die_on_unknown_options();

  if(ufc_hart_fname == NULL)
    {
      n = strlen(output_directory) + 1 + strlen(jobname);
      if(ufc_label != NULL) n += 1 + strlen(ufc_label) + 4;
      str = cart_alloc(char,n);

      if(ufc_label == NULL)
	{
	  sprintf((char *)str,"%s/%s.dh",output_directory,jobname);
	}
      else
	{
	  sprintf((char *)str,"%s/%s_%s.dh",output_directory,jobname,ufc_label);
	}
      ufc_hart_fname = str;
    }
}


void read_file()
{
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
}


void write_file()
{
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
}
