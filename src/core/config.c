#include "config.h"

#include <limits.h>
#include <stdio.h>
#include <string.h>

#include "agn.h"
#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "density.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "io.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "plugin.h"
#include "refinement.h"
#include "rt.h"
#include "starformation.h"
#include "system.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#ifdef _OPENMP
#include <omp.h>
#endif


extern int step;


void config_verify();
void config_dependent_parameters();


void config_create_file(const char *filename)
{
  FILE *f;
  
  if(local_proc_id != MASTER_NODE) return;

  if(filename == NULL)
    {
      f = stdout;
    }
  else
    {
      f = fopen(filename,"w");
      if(f == NULL)
	{
	  cart_error("Unable to open %s for writing!",filename);
	}
    }

  control_parameter_print(f,1);

  if(filename != NULL)
    {
      fclose(f);
    }
}


void config_read_file(const char *filename) {
  FILE *configfile;
  char line[1024];
  char *tag, *value, *p;
  int i; 
  int bad_config = 0;

  configfile = fopen(filename,"r");
  if(configfile == NULL)
    {
      cart_error("Unable to load configuration file %s",filename);
    }

  while(fgets(line,1024,configfile) != NULL)
    {
      /* translate all tabs/carriage return to spaces */
      for(i=0; i<strlen(line); i++)
	{
	  if(line[i]=='\t' || line[i]=='\n')
	    {
	      line[i] = ' ';
	    }
	}

      /* skip over leading spaces */
      p = line;
      while(*p == ' ') p++;

      /* read if this is not a comment line */
      if(p[0] != '#')
	{
	  tag = p;

	  /* Find where it ends */
	  while(*p != ' ') p++;

	  /* null-terminate tag string */
	  *p = 0;
	  p++;

	  /* skip spaces between tag and output */
	  while (*p == ' ') p++;

	  value = p;

	  /* remove trailing spaces */
	  i = strlen(value) - 1;
	  while(i>=0 && value[i] == ' ')
	    {
	      value[i] = 0;
	      i--;
	    }

	  /* read parameters */
	  if(tag[0] != 0) /* ignore empty strings */
	  {
	    if(control_parameter_read(tag,value) != 0)
	      {
		cart_debug("Unrecognized tag [%s] in parameter file %s", tag, filename );
        bad_config = 1;
	      }
	  }
	}
    }

  fclose(configfile);

  if ( bad_config ) {
	cart_error("There were errors parsing parameter file %s", filename );
  }
  
  config_dependent_parameters();

  config_verify();
}


void config_print_to_file(const char *filename, int append)
{
  const char *title_sep = "******************************************\n";
  FILE *f;
  char pathname[256];
  char hostname[256];
  char task_hostnames[MAX_PROCS][256]; 
  int local_pid, pid[MAX_PROCS];
  int proc, n;
  int d_flag[999];

  system_get_host_name( hostname, 256 );
  local_pid = system_get_pid(); 

  MPI_Gather( hostname, 256, MPI_CHAR, task_hostnames, 256, MPI_CHAR, MASTER_NODE, mpi.comm.run );
  MPI_Gather( &local_pid, 1, MPI_INT, pid, 1, MPI_INT, MASTER_NODE, mpi.comm.run );

  if(local_proc_id != MASTER_NODE) return;

  if(filename == NULL)
    {
      f = stdout;
    }
  else
    {
      sprintf(pathname,"%s/%s",logfile_directory,filename);

      f = fopen(pathname,append?"a":"w");
      if(f == NULL)
	{
	  cart_error("Unable to open %s for writing",pathname);
	}
    }

  if(append)
    {
      fprintf(f,"### NEW RUN\n");
      fprintf(f,title_sep);
      fprintf(f,"Number of MPI tasks: %d\n", num_procs );

#ifdef _OPENMP
#pragma omp parallel shared(n)
      n = omp_get_num_threads();
#else
      n = 1;
#endif
      fprintf(f,"Number of OpenMP threads (per-MPI task): %d\n",n);

      for ( proc = 0; proc < num_procs; proc++ ) {
	fprintf(f,"mpi task %3u hostname %s:%d\n", proc, task_hostnames[proc], pid[proc] );
      }

      fprintf(f,"UTC time/date: %s",system_get_time_stamp(1));
      fprintf(f,"Local time/date: %s",system_get_time_stamp(0));

      fprintf(f,"Step: %d, time: %lg",step,tl[min_level]);
#ifdef COSMOLOGY
      fprintf(f,", auni: %lg, abox: %lg",auni[min_level],abox[min_level]);
#endif
      fprintf(f,"\n");
    }

  fprintf(f,"SVN Branch: %s\n", SVNBRANCH );
  fprintf(f,"SVN Revision: %s\n", SVNREVISION );
  fprintf(f,title_sep);

  fprintf(f,"\n");
  fprintf(f,"   GLOBAL SETTINGS\n");
  fprintf(f,title_sep);

#include "../../config/list_defines.h"

  fprintf(f,"Primary settings:\n");
  PRINT(HYDRO);
  PRINT(GRAVITY); 
  PRINT(COOLING);
  PRINT(STAR_FORMATION);
  PRINT(COSMOLOGY);
  PRINT(PARTICLES);
  PRINT(REFINEMENT);
  PRINT(RADIATIVE_TRANSFER);
  PRINT(ELECTRON_ION_NONEQUILIBRIUM);
//  PRINT(TURBULENT_ENERGY);
//  PRINT(COSMIC_RAY_ENERGY);

#ifdef STAR_FORMATION
  fprintf(f,"Star formation settings:\n");
  PRINT(ENRICHMENT);
  PRINT(ENRICHMENT_SNIa);
  PRINT(STAR_PARTICLE_TYPES);
  PRINT(AGN);
#endif /* STAR_FORMATION */

#ifdef RADIATIVE_TRANSFER 
  fprintf(f,"Radiative transfer settings:\n");
  PRINT(RT_TRANSFER);
  PRINT(RT_TRANSFER_METHOD);
  PRINT(RT_UV);
  PRINT(RT_CHEMISTRY);
  PRINT(RT_EXACT_EOS);
  PRINT(RT_HIGH_DENSITY);
  PRINT(RT_XRAYS);
  PRINT(RT_LWBANDS);
  PRINT(RT_LYMAN_ALPHA_HEATING);
  PRINT(RT_PAH_CR);
  PRINT(RT_EXTERNAL_BACKGROUND);
#endif /* RADIATIVE_TRANSFER */

  fprintf(f,"Other settings:\n");
  PRINT_ALL;

  fprintf(f,"\n");
  fprintf(f,"   MESH PARAMETERS\n");
  fprintf(f,title_sep);
#ifdef num_root_grid_refinements
  fprintf(f,"num_root_grid_refinements=\t%-d\n",num_root_grid_refinements);
#endif
#ifdef num_refinement_levels
  fprintf(f,"num_refinement_levels=\t\t%-d\n",num_refinement_levels);
#endif
#ifdef num_octs
  fprintf(f,"num_octs=\t\t\t%-d\n",num_octs);
#endif
#ifdef num_particles
  fprintf(f,"num_particles=\t\t\t%-d\n",num_particles);
#endif
#ifdef num_star_particles
  fprintf(f,"num_star_particles=\t\t%-d\n",num_star_particles);
#endif

  fprintf(f,"\n");
  fprintf(f,"   CONTROL PARAMETERS\n");
  fprintf(f,title_sep);

  control_parameter_print(f,0);

  fprintf(f,"\n");
  fprintf(f,"!!!\n");
  fprintf(f,"!!! CHANGING THESE PARAMETERS WILL BREAK THE CODE\n");
  fprintf(f,"!!!  THEY ARE LISTED FOR DIAGNOSTIC PURPOSE ONLY\n");
  fprintf(f,"!!!\n");

  control_parameter_print_hidden(f,0);

  if(filename != NULL)
    {
      fclose(f);
    }
}


void config_append_units_to_file(const char *filename)
{
  const char *title_sep = "******************************************\n";
  FILE *f;
  char pathname[256];
  
  if(local_proc_id != MASTER_NODE) return;

  if(filename == NULL)
    {
      f = stdout;
    }
  else
    {
      sprintf(pathname,"%s/%s",logfile_directory,filename);
      f = fopen(pathname,"a");
      if(f == NULL)
	{
	  cart_error("Unable to open %s for writing!",pathname);
	}
    }

  fprintf(f,"\n");
  fprintf(f,"   UNITS\n");
  fprintf(f,title_sep);

#ifdef COSMOLOGY

  fprintf(f,"Cosmological parameter:\n");
  fprintf(f,"   H0:      %lg km/s/Mpc\n",cosmology->h*100);
  fprintf(f,"   OmegaM:  %lg\n",cosmology->OmegaM);
  fprintf(f,"   OmegaD:  %lg\n",cosmology->OmegaD);
  fprintf(f,"   OmegaB:  %lg\n",cosmology->OmegaB);
  fprintf(f,"   OmegaL:  %lg\n",cosmology->OmegaL);
  fprintf(f,"   OmegaK:  %lg\n",cosmology->OmegaK);
  fprintf(f,"   DeltaDC: %lg\n",cosmology->DeltaDC);
  if(cosmology->flat)  fprintf(f,"   This cosmology is flat.\n");
  fprintf(f,"\n");
  fprintf(f,"   Box size:\t \t\t%-g CHIMP\n",box_size);
  fprintf(f,"\n");

#endif /* COSMOLOGY */

  fprintf(f,"Primary units:\n");
  fprintf(f,"   Mass:      %lg g\n",primary_units->mass);
  fprintf(f,"   Time:      %lg s\n",primary_units->time);
  fprintf(f,"   Length:    %lg cm\n",primary_units->length);

  if(filename != NULL)
    {
      fclose(f);
    }
}


void config_allocate_data(float memory_fraction_mesh);

void config_init()
{
  config_init_plugins();

  config_init_io();
  config_init_load_balance();
  config_init_times();

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
  config_init_density();
#endif

#ifdef HYDRO
  config_init_hydro();
#endif

#ifdef GRAVITY
  config_init_gravity();
#endif

#ifdef REFINEMENT
  config_init_refinement();
#endif

#ifdef STAR_FORMATION
  config_init_star_formation();
#ifdef AGN
  config_init_agn();
#endif /* AGN */
#endif

#ifdef RADIATIVE_TRANSFER
  rtConfigInit();
#endif /* RADIATIVE_TRANSFER */

  config_init_parallel();
}

void config_dependent_parameters()
{
#ifdef STAR_FORMATION
  config_dependent_star_formation_parameters();
#endif
}

void config_verify()
{
  config_verify_plugins();

  config_verify_io();
  config_verify_load_balance();
  config_verify_times();

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
  config_verify_density();
#endif

#ifdef HYDRO
  config_verify_hydro();
#endif

#ifdef REFINEMENT
  config_verify_refinement();
#endif

#ifdef STAR_FORMATION
  config_verify_star_formation();
#ifdef AGN
  config_verify_agn();
#endif /* AGN */  
#endif

#ifdef RADIATIVE_TRANSFER
  rtConfigVerify();
#endif /* RADIATIVE_TRANSFER */

  config_verify_parallel();
}

