#include "config.h"


#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "cosmology.h"
#include "times.h"
#include "units.h"

#include "../core/cell_buffer.h"
#include "../core/density.h"
#include "../core/hydro.h"
#include "../core/io.h"
#include "../core/iterators.h"
#include "../core/load_balance.h"
#include "../core/parallel.h"
#include "../core/particle.h"
#include "../core/plugin.h"
#include "../core/rt.h"
#include "../core/tree.h"

#include "start_analysis.h"

extern int max_first_species_id;


void config_append_units_to_file(const char *filename);
void config_print_to_file(const char *filename, int append);


/*
//  Initialize what is needed
*/
int  init_da_set_ = 0;
void init_da_exe_()
{
  int level;

  if(!init_da_set_)
    {
      init_da_set_ = 1;

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
      /*
      //  This trick ensures that cell_first_species_mass(c) is now stellar density
      */
      max_first_species_id = 0;
      for(level=min_level; level<=max_level; level++)
	{
	  cart_debug("assigning density on level %u", level );
	  assign_density( level );
	}
      max_first_species_id = -1; /* This will automatically reset it at the next density assignment call. */
    }
#endif /* GRAVITY || RADIATIVE_TRANSFER */
}


int  init_td_set_ = 0;
void init_td_exe_()
{
  MESH_RUN_DECLARE(level,cell);

  if(!init_td_set_)
    {
      init_td_set_ = 1;

      init_da_exe_();

#ifdef GRAVITY
      MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,cell_volume_inverse)
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
      
      cell_total_density(cell) = cell_total_mass(cell)*cell_volume_inverse[level] + 1.0;

      MESH_RUN_OVER_CELLS_OF_LEVEL_END;
      MESH_RUN_OVER_LEVELS_END;

      cart_debug("Total density is assigned.");
#endif /* GRAVITY */
    }
}


int  init_sd_set_ = 0;
void init_sd_exe_()
{
  MESH_RUN_DECLARE(level,cell);

  if(!init_sd_set_)
    {
      init_sd_set_ = 1;

      init_da_exe_();

#if defined(GRAVITY) && defined(STARFORM)
      MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,cell_volume_inverse)
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
      
      cell_stellar_density(cell) *= cell_volume_inverse[level];

      MESH_RUN_OVER_CELLS_OF_LEVEL_END;
      MESH_RUN_OVER_LEVELS_END;

      cart_debug("Stellar density is assigned.");
#endif /* GRAVITY && STARFORM */
    }
}


int  init_rt_set_ = 0;
void init_rt_exe_()
{
#ifdef RADIATIVE_TRANSFER
  if(!init_rt_set_)
    {
      init_rt_set_ = 1;

      rtInitStep(0.0);
      init_da_exe_();

      cart_debug("RT is initialized.");
    }
#endif
}


void init_reset_()
{
  init_da_set_ = 0;
  init_td_set_ = 0;
  init_sd_set_ = 0;
  init_rt_set_ = 0;
}


void init_all_()
{
  init_td_exe_();
  init_sd_exe_();
  init_rt_exe_();
}

struct init_t init_internal_ = { init_all_, init_td_exe_, init_sd_exe_, init_rt_exe_ };
const struct init_t *init = &init_internal_;


struct snapshot_t snapshot_internal_ = { 0, 0 };
const struct snapshot_t *snapshot = &snapshot_internal_;


extern int current_step_level;

/*
//  Default plugin callback
*/
void run_output()
{
}


/*
//  Analysis engine
*/
int main_analysis(int argc, const char *argv[]);

void run(int restart, const char *restart_label)
{
  int i, ret, level, full_mode;
  double aexp;
  char str[999];
  int argc, num_snapshots;
  const char **argv;

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Entering analysis mode...");
    }

  if(num_options>0 && (strcmp(options[0],"-f")==0 || strcmp(options[0],"--fast")==0))
    {
      cart_debug("Starting in a fast mode, without initializations...");
      num_options--;
      options++;
      full_mode = 0;
    }
  else
    {
      full_mode = 1;
    }

  for(i=0; i<num_options; i++)
    {
      if(options[i][0] == '-') break;
    }
  num_snapshots = i;

  if(num_snapshots == 0)
    {
      cart_debug("At least one snapshot needs to be specified after all ART options:");
      cart_error("analysis [ART options] <aexp1> [aexp2 ...] [analysis options]");
    }

  argc = num_options - num_snapshots;
  argv = options + num_snapshots;

  snapshot_internal_.Number = num_snapshots;

  for(i=0; i<num_snapshots; i++)
    {
#ifdef COSMOLOGY
      aexp = -1.0;
      if(sscanf(options[i],"%lf%c",&aexp,str)!=1 || aexp<1.0e-10 || aexp>1.1)
	{
	  cart_debug("Invalid scale factor value %s; skipping snapshot #%d",options[i],i);
	  continue;
	}
#endif /* COSMOLOGY */
      
#ifdef COSMOLOGY
      strcpy(str+1,options[i]);
      str[0] = 'a';
#else
      strcpy(str,options[i]);
#endif /* COSMOLOGY */

      snapshot_internal_.Current = i;

      if(num_snapshots > 1)
	{
#ifdef COSMOLOGY
	  cart_debug("Analysing snapshot #%d at aexp=%6.4f",i+1,aexp);
#else
	  cart_debug("Analysing snapshot #%d named %s",i+1,str);
#endif /* COSMOLOGY */
	  current_step_level = 3;
	}
      else
	{
	  current_step_level = -1;
	}

      if(buffer_enabled) destroy_cell_buffer();
#ifdef PARTICLES
      init_particles();
#endif /* PARTICLES */

      read_restart(str);
      load_balance();

#ifdef COSMOLOGY
      abox[min_level] = abox_from_tcode(tl[min_level]);
      auni[min_level] = auni_from_tcode(tl[min_level]);
#endif /* COSMOLOGY */

      for(level=min_level+1; level<=max_level; level++)
	{
	  tl[level] = tl[min_level];
#ifdef COSMOLOGY
	  abox[level] = abox[min_level];
	  auni[level] = auni[min_level];
#endif /* COSMOLOGY */
	}

      init_reset_();
      if(full_mode) init->All();

      config_append_units_to_file("config.log");
      config_print_to_file("history.log",1);

      ret = main_analysis(argc,argv);

      current_step_level = -1;

      if(ret != 0)
	{
	  cart_debug("main_analysis returned with a non-zero return code %d. Aborting...",ret);
	  break;
	}
    }
}


