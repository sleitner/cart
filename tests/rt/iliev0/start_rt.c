#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "defs.h"
#include "tree.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "viewdump.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "gravity.h"
#include "density.h"
#include "io.h"
#include "auxiliary.h"

#include "rt_solver.h"
#include "rt_utilities.h"


#define N50             0.0
#define T_i             1.0e2


void rt_initial_conditions( int cell )
{
  cell_gas_density(cell) = 1.0;
  cell_momentum(cell,0) = 0.0;
  cell_momentum(cell,1) = 0.0;
  cell_momentum(cell,2) = 0.0;
  cell_gas_gamma(cell) = (5.0/3.0);

  cell_gas_internal_energy(cell) = 1.5*T_i*wmu/T0*(rt_XH+rt_XHe);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)/1.5;
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);

  cell_HI_density(cell) = rt_XH;
  cell_HII_density(cell) = 0.0;
  cell_HeI_density(cell) = rt_XHe;
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
      fprintf(F,"%e %e %e %e\n",(float)(t0*tl[min_level]),cell_HI_density(0)/cell_gas_density(0),rtTemInK(0),cell_H2_density(0)/cell_gas_density(0));
      fclose(F);
    }
}


void init_run() {
	int i;
	int level;
	int ioct;
	int num_level_cells;
	int *level_cells;
        int all_hydro_vars[num_hydro_vars];
	char filename[128];
	int min_index, max_index;

	/* set units */
	Lbox = 2*15e-3;
	hubble = 1.0;
	Omega0 = 1.0e-0/(1.123e-5*hubble*hubble);
	Omegab0 = Omega0;
	OmegaL0 = 0.0;

        /* create array with all hydro variable indices */
        for ( i = 0; i < num_hydro_vars; i++ ) {
                all_hydro_vars[i] = HVAR_GAS_DENSITY + i;
        }

	for ( i = 0; i < nDim; i++ ) {
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

	for ( level = min_level; level <= max_level; level++ ) {
		cart_debug("updating level %u", level );
		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	}

	cart_debug("done updating initial conditions");

	/* set time variables */
	tl[min_level] = 0.0;

	dtl[min_level] = 1.0e-5/t0;
	/*choose_timestep( &dtl[min_level] );*/

#ifdef COSMOLOGY
	aexp[min_level] = b2a( tl[min_level] );
#else
	aexp[min_level] = 1.0;
#endif

	for ( level = min_level+1; level <= max_level; level++ ) {
		dtl[level] = 0.5*dtl[level-1];
		tl[level] = tl[min_level];
		aexp[level] = aexp[min_level];		
	}

	cart_debug("done with initialization");

	check_map();

	run_output();
}
