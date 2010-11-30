
#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "refinement.h"
#include "rt_solver.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#ifdef RT_DEBUG
#include "rt_debug.h"
#endif

const float T_i = 1.0e2;


void units_set_art(double OmegaM, double h, double Lbox);


void rt_initial_conditions( int cell )
{
  cell_gas_density(cell) = 1.0;
  cell_momentum(cell,0) = 0.0;
  cell_momentum(cell,1) = 0.0;
  cell_momentum(cell,2) = 0.0;
  cell_gas_gamma(cell) = (5.0/3.0);

  cell_gas_internal_energy(cell) = T_i/units->temperature/(constants->gamma-1);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(constants->gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);

  cell_HI_density(cell) = 1.0;
  cell_HII_density(cell) = 0.0;
  cell_HeI_density(cell) = 1.0e-10;
  cell_HeII_density(cell) = 0.0;
  cell_HeIII_density(cell) = 0.0;
  cell_H2_density(cell) = 0.0;
}


void set_rt_initial_conditions(void)
{
  MESH_RUN_DECLARE(level,cell);

  _MaxLevel = min_level;

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  rt_initial_conditions(cell);
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  for(level=_MaxLevel; level>=min_level; level--)
    {
      hydro_split_update(level);
    }
}


void run_output()
{
  FILE *F;

  if(local_proc_id == MASTER_NODE)
    {
      if(step == 0) F = fopen("evol.res","w"); else F = fopen("evol.res","a"); 
      fprintf(F,"%e %e %e %e\n",(float)(units->time*tl[min_level]/constants->yr),cell_HI_density(0)/cell_gas_density(0),units->temperature*rtTem(0),cell_H2_density(0)/cell_gas_density(0));
      fclose(F);
    }
}


void init_run()
{
  int i;
  int level;
  float hubble;

  /* set units */
  hubble = 1;
  units_set_art(1.0e-0/(1.123e-5*hubble*hubble),hubble,2*15e-3/hubble);

  units_reset();
  units_update(min_level);

  for ( i = 0; i < nDim; i++ )
    {
      refinement_volume_min[i] = 0.0;
      refinement_volume_max[i] = (double)num_grid;
    }

  cart_debug("in init");

  /* build buffer */
  build_cell_buffer();
  cart_debug("built cell buffer");
  repair_neighbors();

  check_map();

  cart_debug("repaired neighbors");

  cart_debug("setting initial conditions");

  set_rt_initial_conditions();

  cart_debug("set initial conditions");

  for ( level = min_level; level <= max_level; level++ )
    {
      cart_debug("updating level %u", level );
      update_buffer_level( level, all_hydro_vars, num_hydro_vars );
    }

  cart_debug("done updating initial conditions");

  /* set time variables */
  tl[min_level] = 0.0;

  dtl[min_level] = 1.0e-5*constants->yr/units->time;
  /*choose_timestep( &dtl[min_level] );*/

  for ( level = min_level+1; level <= max_level; level++ )
    {
      dtl[level] = 0.5*dtl[level-1];
      tl[level] = tl[min_level];
    }

  cart_debug("done with initialization");

  check_map();

  run_output();
  /*
  //  Debugging parameters
  */
#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = 1;
  rt_debug.Stop = 0;
  cell_center_position(1,rt_debug.Pos);
#endif
#endif
}
