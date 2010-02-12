#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "density.h"
#include "hydro.h"
#include "io.h"
#include "load_balance.h"
#include "parallel.h"
#include "refinement.h"
#include "rt_solver.h"
#include "starformation.h"
#include "system.h"
#include "timestep.h"
#include "units.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void config_verify();


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

  config_verify();
}


void config_print_to_file(const char *filename, int restart)
{
  const char *title_sep = "******************************************\n";
  FILE *f;
  char pathname[256];
  char hostname[256];
  char task_hostnames[MAX_PROCS][256]; 
  int local_pid, pid[MAX_PROCS];
  int proc;
  int n;

  system_get_host_name( hostname, 256 );
  local_pid = system_get_pid(); 

  MPI_Gather( hostname, 256, MPI_CHAR, task_hostnames, 256, MPI_CHAR, MASTER_NODE, MPI_COMM_WORLD );
  MPI_Gather( &local_pid, 1, MPI_INT, pid, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

  if(local_proc_id != MASTER_NODE) return;

  if(filename == NULL)
    {
      f = stdout;
    }
  else
    {
      sprintf(pathname,"%s/%s",logfile_directory,filename);

      if ( restart ) {
	      f = fopen(pathname,"w");
      } else {
          f = fopen(pathname,"a");
      }
      if(f == NULL)
	{
	  cart_error("Unable to open %s for writing!",pathname);
	}
    }

  fprintf(f,"SVN Branch: %s\n", SVNBRANCH );
  fprintf(f,"SVN Revision: %s\n", SVNREVISION );
  fprintf(f,title_sep);

  fprintf(f,"\n");
  fprintf(f,"   GLOBAL SETTINGS\n");
  fprintf(f,title_sep);

  fprintf(f,"Primary settings:\n");
#ifdef HYDRO
  fprintf(f,"   HYDRO\n");
#endif
#ifdef GRAVITY 
  fprintf(f,"   GRAVITY\n");
#endif
#ifdef COOLING 
  fprintf(f,"   COOLING\n");
#endif
#ifdef STARFORM 
  fprintf(f,"   STARFORM\n");
#endif
#ifdef COSMOLOGY 
  fprintf(f,"   COSMOLOGY\n");
#endif
#ifdef PARTICLES 
  fprintf(f,"   PARTICLES\n");
#endif
#ifdef REFINEMENT 
  fprintf(f,"   REFINEMENT\n");
#endif
#ifdef RADIATIVE_TRANSFER 
  fprintf(f,"   RADIATIVE_TRANSFER\n");
#else
#ifdef NO_METALCOOLING 
  fprintf(f,"   NO_METALCOOLING\n");
#endif
#endif

#ifdef STARFORM
  fprintf(f,"Star formation settings:\n");
#ifdef ENRICH 
  fprintf(f,"   ENRICH\n");
#endif
#ifdef ENRICH_SNIa
  fprintf(f,"   ENRICH_SNIa\n");
#endif
#ifdef FEEDBACK 
  fprintf(f,"   FEEDBACK\n");
#endif
#ifdef FEEDBACK_SNIa
  fprintf(f,"   FEEDBACK_SNIa\n");
#endif
#ifdef STELLARMASSLOSS 
  fprintf(f,"   STELLARMASSLOSS\n");
#endif
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER 
  fprintf(f,"Radiative transfer settings:\n");
#ifdef RT_TABLES 
  fprintf(f,"   RT_TABLES\n");
#endif
#ifdef RT_TRANSFER 
  fprintf(f,"   RT_TRANSFER\n"); 
#endif
#ifdef RT_TRANSFER_METHOD 
  fprintf(f,"   RT_TRANSFER_METHOD=%-d\n",RT_TRANSFER_METHOD);
#endif
#ifdef RT_CHEMISTRY 
  fprintf(f,"   RT_CHEMISTRY\n"); 
#endif 
#ifdef RT_MONOATOMIC 
  fprintf(f,"   RT_MONOATOMIC\n"); 
#endif
#ifdef RT_HIGH_DENSITY 
  fprintf(f,"   RT_HIGH_DENSITY\n"); 
#endif
#ifdef RT_XRAYS 
  fprintf(f,"   RT_XRAYS\n"); 
#endif
#ifdef RT_DUST 
  fprintf(f,"   RT_DUST\n"); 
#endif
#ifdef RT_LWBANDS 
  fprintf(f,"   RT_LWBANDS\n"); 
#endif
#ifdef RT_LYMAN_ALPHA_HEATING 
  fprintf(f,"   RT_LYMAN_ALPHA_HEATING\n"); 
#endif
#ifdef RT_PAH_CR 
  fprintf(f,"   RT_PAH_CR\n"); 
#endif
#ifdef RT_SIGNALSPEED_TO_C 
  fprintf(f,"   RT_SIGNALSPEED_TO_C\n"); 
#endif
#ifdef RT_EXTERNAL_BACKGROUND 
  fprintf(f,"   RT_EXTERNAL_BACKGROUND=%-d\n",RT_EXTERNAL_BACKGROUND);
#endif
#ifdef RT_CFI 
  fprintf(f,"   RT_CFI=%-d\n",RT_CFI);
#endif
#ifdef RT_DEBUG 
  fprintf(f,"   RT_DEBUG\n"); 
#endif
#ifdef RT_OUTPUT 
  fprintf(f,"   RT_OUTPUT\n"); 
#endif
#ifdef RT_PARALLEL_NUM_OPENMP_BUFFERS 
  fprintf(f,"   RT_PARALLEL_NUM_OPENMP_BUFFERS=%-d\n",RT_PARALLEL_NUM_OPENMP_BUFFERS);
#endif
#ifdef RT_PARALLEL_USE_MPI 
  fprintf(f,"   RT_PARALLEL_USE_MPI\n"); 
#endif
#ifdef RT_INTERPOLLOG 
  fprintf(f,"   RT_INTERPOLLOG\n"); 
#endif
#ifdef RT_NARROWTABLE 
  fprintf(f,"   RT_NARROWTABLE\n"); 
#endif
#ifdef RT_8SPECIES 
  fprintf(f,"   RT_8SPECIES \n");
#endif
#ifdef RT_TRANSFER_FLUX_CONSERVING 
  fprintf(f,"   RT_TRANSFER_FLUX_CONSERVING\n"); 
#endif
#ifdef RT_H2_RATE
  fprintf(f,"   RT_H2_RATE=%-d\n",RT_H2_RATE);
#endif
#ifdef RT_H2_COOL
  fprintf(f,"   RT_H2_COOL=%-d\n",RT_H2_COOL);
#endif
#ifdef RT_VARIABLE_PRATES 
  fprintf(f,"   RT_VARIABLE_PRATES\n"); 
#endif
#ifdef RT_VARIABLE_RFIELD 
  fprintf(f,"   RT_VARIABLE_RFIELD\n");
#endif
#endif /* RADIATIVE_TRANSFER */

  fprintf(f,"Secondary settings:\n");
#ifdef LAPIDUS 
  fprintf(f,"   LAPIDUS\n");
#endif
#ifdef DENSGRADSMOOTH 
  fprintf(f,"   DENSGRADSMOOTH\n");
#endif
#ifdef PRESSURE_FLOOR 
  fprintf(f,"   PRESSURE_FLOOR\n");
#endif
#ifdef ADVECT_SPECIES
  fprintf(f,"   ADVECT_SPECIES\n");
#endif
#ifdef PRESSURELESS_FLUID
  fprintf(f,"   PRESSURELESS_FLUID\n");
#endif
#ifdef ELECTRON_ION_NONEQUILIBRIUM
  fprintf(f,"   ELECTRON_ION_NONEQUILIBRIUM\n");
#endif
#ifdef DEBUG 
  fprintf(f,"   DEBUG\n");
#endif
#ifdef DEBUG_MEMORY_USE
  fprintf(f,"   DEBUG_MEMORY_USE\n");
#endif
#ifdef OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION 
  fprintf(f,"   OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION\n"); 
#endif
#ifdef OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID 
  fprintf(f,"   OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID\n");
#endif

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

  fprintf(f,"\n");
  fprintf(f,"   SYSTEM SETTINGS\n");
  fprintf(f,title_sep);

  fprintf(f,"\n");
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
  fprintf(f,"   Mass:      %lg g\n",primary_units->mass/cgs->g);
  fprintf(f,"   Time:      %lg s\n",primary_units->time/cgs->s);
  fprintf(f,"   Length:    %lg cm\n",primary_units->length/cgs->cm);

#ifdef LEGACY_UNITS
  fprintf(f,"\n");
  fprintf(f,"Legacy units:\n");
  fprintf(f,"   H0 = %le\n",legacy_units->H0);
  fprintf(f,"   r0 = %le\n",legacy_units->r0);
  fprintf(f,"   t0 = %le\n",legacy_units->t0);
  fprintf(f,"   v0 = %le\n",legacy_units->v0);
  fprintf(f,"   rho0 = %le\n",legacy_units->rho0);
  fprintf(f,"   den0 = %le\n",legacy_units->den0);
  fprintf(f,"   P0 = %le\n",legacy_units->P0);
  fprintf(f,"   T0 = %le\n",legacy_units->T0);
  fprintf(f,"   E0 = %le\n",legacy_units->E0);
  fprintf(f,"   M0 = %le\n",legacy_units->M0);
#endif /* LEGACY_UNITS */

  if(filename != NULL)
    {
      fclose(f);
    }
}


void config_init()
{
  /* This MUST be the first call */
  config_init_units();

  config_init_io();
  config_init_load_balance();
  config_init_timestep();

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
  config_init_density();
#endif

#ifdef HYDRO
  config_init_hydro();
#endif

#ifdef REFINEMENT
  config_init_refinement();
#endif

#ifdef STARFORM
  config_init_star_formation();
#endif

#ifdef RADIATIVE_TRANSFER
  rtConfigInit();
#endif /* RADIATIVE_TRANSFER */

  config_init_parallel();
}


void config_verify()
{
  /* This MUST be the first call */
  config_verify_units();

  config_verify_io();
  config_verify_load_balance();
  config_verify_timestep();

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
  config_verify_density();
#endif

#ifdef HYDRO
  config_verify_hydro();
#endif

#ifdef REFINEMENT
  config_verify_refinement();
#endif

#ifdef STARFORM
  config_verify_star_formation();
#endif

#ifdef RADIATIVE_TRANSFER
  rtConfigVerify();
#endif /* RADIATIVE_TRANSFER */

  config_verify_parallel();
}

