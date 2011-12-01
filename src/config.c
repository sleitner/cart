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

#include "../config/list_defines.h"

  fprintf(f,"Primary settings:\n");
  PRINT(HYDRO);
  PRINT(GRAVITY); 
  PRINT(COOLING);
  PRINT(STARFORM);
  PRINT(COSMOLOGY);
  PRINT(PARTICLES);
  PRINT(REFINEMENT);
  PRINT(RADIATIVE_TRANSFER);
  PRINT(ELECTRON_ION_NONEQUILIBRIUM);

#ifdef STARFORM
  fprintf(f,"Star formation settings:\n");
  PRINT(ENRICH);
  PRINT(ENRICH_SNIa);
  PRINT(STAR_PARTICLE_TYPES);
  PRINT(AGN);
#endif /* STARFORM */

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
  fprintf(f,"   Mass:      %lg g\n",primary_units->mass/constants->g);
  fprintf(f,"   Time:      %lg s\n",primary_units->time/constants->s);
  fprintf(f,"   Length:    %lg cm\n",primary_units->length/constants->cm);

#ifdef CHECK_LEGACY_UNITS
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
#endif /* CHECK_LEGACY_UNITS */

  if(filename != NULL)
    {
      fclose(f);
    }
}


void config_allocate_data(float memory_fraction_mesh);

void config_init()
{
#ifndef STATIC_MESH_DATA
  config_allocate_data(0.6);
#endif /* STATIC_MESH_DATA */

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

#ifdef GRAVITY
  config_init_gravity();
#endif

#ifdef REFINEMENT
  config_init_refinement();
#endif

#ifdef STARFORM
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
#ifdef AGN
  config_verify_agn();
#endif /* AGN */  
#endif

#ifdef RADIATIVE_TRANSFER
  rtConfigVerify();
#endif /* RADIATIVE_TRANSFER */

  config_verify_parallel();
}


#ifndef STATIC_MESH_DATA
extern int OCT_ARRAY(cells_to_refine);
extern const int is_running;

double config_allocate_mesh_helper(int mode, int nocts);

/*
//  Allocate data if we are running in a dynamic mode
*/
void config_allocate_data(float memory_fraction_mesh)
{
  const double GB = 1024.0*1024.0*1024.0;
  double maxmem, taskmem, mem_per_oct, mem;
  const char *str;
  int min_octs = (num_root_cells >> nDim)/num_procs;
  size_t max_octs;
  int do_not_zero_data;

  cart_assert(cell_vars == NULL);
  cart_assert(num_octs == 0);

  /*
  //  Find max available (or requested) memory per task
  */
  taskmem = system_get_available_memory()/tasks_per_node;

  str = extract_option1("max-memory","m",NULL);
  if(str != NULL)
    {
      if(sscanf(str,"%lg",&maxmem) != 1)
        {
          cart_error("A value for the max memory per task (in GB) must be specified.");
        }

      if(maxmem < 1.0e-3)
        {
          cart_error("A value of %lg for the max memory per task is not valid.",maxmem);
        }
    }
  else maxmem = 1.0e30;

  maxmem *= GB;
  if(maxmem > taskmem) maxmem = taskmem;

  /*
  //  Specify how many octs we want.
  */
  str = extract_option1("num-octs","no",NULL);
  if(str == NULL)
    {
      /*
      //  Use max available memory.
      */
      num_octs = 0;
    }
  else
    {
      if(sscanf(str,"%d",&num_octs) != 1)
        {
          cart_error("A value for the number of octs MUST be specified.");
        }

      if(num_octs!=0 && num_octs<min_octs)
        {
          cart_error("There must be at least %d octs.",min_octs);
        }
    }

  mem_per_oct = config_allocate_mesh_helper(0,1);

  cart_debug("Memory per oct: %lg bytes",mem_per_oct);

  max_octs = (size_t)(0.5+memory_fraction_mesh*maxmem/mem_per_oct);

  /*
  //  Check that we do not wrap an int
  */
  if(max_octs > (INT_MAX >> nDim))
    {
      cart_debug("WARNING: Maximum available number of octs would result in wrapping of an int!!!");
      max_octs = (INT_MAX >> nDim) - 100;
    }

  if(num_octs == 0) num_octs = max_octs;

  if(num_octs > max_octs)
    {
      cart_error("The requested number of octs (%d) is too large; only %d octs are available.",num_octs,max_octs);
    }

  num_cells = (num_octs << nDim);

  cart_debug("Allocating mesh: num_octs=%d, num_cells=%d",num_octs,num_cells);


  do_not_zero_data = is_option_present("fast","f",0);
  mem = config_allocate_mesh_helper(do_not_zero_data?1:2,num_octs);

  cart_debug("Total allocated memory for global data: %5.2lfGB (out of %5.2lfGB)",mem/GB,taskmem/GB);
}


#define ALLOC(array,type,dim,dim2) \
  mem += sizeof(type)*(size_t)dim*dim2; \
  if(mode > 0) array = cart_alloc(type,(size_t)dim*dim2);	\
  if(mode > 1) memset(array,0,sizeof(type)*(size_t)dim*dim2)

double config_allocate_mesh_helper(int mode, int nocts)
{
  int i;
  int *tmpi;
  double *tmpd;
  double mem = 0.0;
  int ncells;

  cart_assert(nocts > 0);
  ncells = (nocts << nDim);

  ALLOC(cell_vars,float,ncells,num_vars);

#ifdef HYDRO
  if(is_running)
    {
      mem += sizeof(float)*ncells*(num_hydro_vars-2+1);
    }
#endif

  ALLOC(cell_child_oct,int,ncells,1);

  ALLOC(oct_parent_cell,int,nocts,1);

  ALLOC(oct_level,int,nocts,1);

  ALLOC(oct_parent_root_sfc,int,nocts,1);

  ALLOC(oct_next,int,nocts,1);

  ALLOC(oct_prev,int,nocts,1);

  ALLOC(oct_neighbors,int*,nocts,1);
  ALLOC(tmpi,int,nocts,num_neighbors);
  if(mode > 0)
    {
      for(i=0; i<nocts; i++)
	{
	  oct_neighbors[i] = tmpi + i*num_neighbors;
	}
    }

  ALLOC(oct_pos,double*,nocts,1);
  ALLOC(tmpd,double,nocts,nDim);
  if(mode > 0)
    {
      for(i=0; i<nocts; i++)
	{
	  oct_pos[i] = tmpd + i*nDim;
	}
    }

#ifdef PARTICLES
  ALLOC(cell_particle_list,int,ncells,1);
#endif /* PARTICLES */

#ifdef REFINEMENT
  ALLOC(cells_to_refine,int,ncells,1);
#endif /* REFINEMENT */

#if defined(HYDRO) && defined(HYDRO_TRACERS)
  ALLOC(cell_tracer_list,int,ncells,1);
#endif /* HYDRO && HYDRO_TRACERS */

  return mem;
}


#endif /* STATIC_MESH_DATA */

