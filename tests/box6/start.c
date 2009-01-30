#include "defs.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include "auxiliary.h"
#include "tree.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "gravity.h"
#include "density.h"
#include "starformation.h"
#include "io.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#include "rt_utilities.h"
#ifdef RT_DEBUG
#include "rt_debug.h"
#endif
#endif

#include "extra/ifrit.h"


void FindMaxVar(int var, float *val, double *pos)
{
  MESH_RUN_DECLARE(level,cell);
  float vMax = -1.0e35;
  int cellMax = 0;

  cart_assert(var>=0 && var<num_vars);

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell) && vMax<cell_var(cell,var))
    {
      vMax = cell_var(cell,var);
      cellMax = cell;
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  *val = vMax;
  if(pos != NULL) cell_position_double(cellMax,pos);
}


void run_output()
{
  const int nbin1 = 256;
#ifdef RADIATIVE_TRANSFER
  int varid[] = { EXT_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, EXT_GAS_TEMPERATURE, EXT_FRACTION+RT_HVAR_OFFSET+5, EXT_CELL_LEVEL, EXT_LOCAL_PROC };
#else
  int varid[] = { HVAR_GAS_DENSITY, HVAR_PRESSURE };
#endif
  int nbin[] = { nbin1, nbin1, nbin1 };
  int nvars = sizeof(varid)/sizeof(int);
  double bb[6], pos[3], dbb;
  float dmax;
  char filename[99];

  bb[0] = 0.0;
  bb[1] = num_grid;
  bb[2] = 0.0;
  bb[3] = num_grid;
  bb[4] = 0.0;
  bb[5] = num_grid;

  sprintf(filename,"%s/out-box.%04d.bin",output_directory,(int)(aexp[0]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

  FindMaxVar(HVAR_GAS_DENSITY,&dmax,pos);

  dbb = nbin1*pow(0.5,(double)max_level);

  bb[0] = pos[0] - 0.5*dbb;
  bb[1] = pos[0] + 0.5*dbb;
  bb[2] = pos[1] - 0.5*dbb;
  bb[3] = pos[1] + 0.5*dbb;
  bb[4] = pos[2] - 0.5*dbb;
  bb[5] = pos[2] + 0.5*dbb;

  sprintf(filename,"%s/out-zoom.%04d.bin",output_directory,(int)(aexp[0]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);
}

void init_run()
{
  int i,j;
  int level;
  int total_cells_per_level[max_level-min_level+1];
  float refmin[nDim];
  float refmax[nDim];
  char filename[256], filename2[256];

#ifdef PARTICLES
  sprintf( filename, "../IC/PMcrd.DAT", output_directory );
  sprintf( filename2, "../IC/PMcrs0.DAT", output_directory );

  restart_load_balance( NULL, filename, filename2 );

  read_particles( filename, filename2, NULL, NULL, 0, NULL );
  cart_debug("read in particles");
#endif

#ifdef HYDRO
  sprintf( filename, "../IC/tr_ic.dat", output_directory );
  read_gas_ic(filename);
  cart_debug("read in gas");

  init_units();

#ifdef RADIATIVE_TRANSFER
  float xH = 1.0 - Y_p;
  float xHe = 0.25*Y_p;
  float xInit[6];

  xInit[1] = 1.2e-5*sqrt(Omega0)/(Omegab0*hubble);
  xInit[0] = xH - xInit[1];
  xInit[4] = xInit[3] = 1.0e-10;
  xInit[2] = xHe - xInit[3] - xInit[4];
  xInit[5] = (aexp[min_level] < 0.0125) ? 7.0e-7/xH : 2.0e-6/xH;

  for(i=0; i<num_cells_per_level[min_level]; i++)
    {
      for(j=0; j<6; j++) cell_vars[i][RT_HVAR_OFFSET+j] = xInit[j]*cell_gas_density(i);
    }
#endif

  hydro_magic( min_level );
  hydro_eos( min_level );
#endif /* HYDRO */

  cart_debug("tl[min_level] = %f", tl[min_level] );
  cart_debug("aexp[min_level] = %f", aexp[min_level] );

  dtl[min_level] = 0.0;
  choose_timestep( &dtl[min_level] );

#ifdef PARTICLES
  for ( i = 0; i < nDim; i++ ) {
    refmin[i] = num_grid+1.0;
    refmax[i] = -1.0;

    for ( j = 0; j < num_particles; j++ )
      {
	if ( particle_level[j] != FREE_PARTICLE_LEVEL && particle_id[j] < particle_species_indices[1] )
	  {
	    if ( particle_x[j][i] < refmin[i] )
	      {
		refmin[i] = particle_x[j][i];
	      }

	    if ( particle_x[j][i] > refmax[i] )
	      {
		refmax[i] = particle_x[j][i];
	      }
	  }
      }
  }
#else
  for ( i = 0; i < nDim; i++ )
    {
      refmin[i] = -1.0;
      refmax[i] = num_grid+1;
    }
#endif /* PARTICLES */

  MPI_Allreduce( refmin, refinement_volume_min, nDim, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD );
  MPI_Allreduce( refmax, refinement_volume_max, nDim, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD );

  for ( i = 0; i < nDim; i++ )
    {
      cart_debug("refinement_volume[%u] = %e %e", i, refinement_volume_min[i], refinement_volume_max[i] );
    }

#ifdef STARFORM
  for ( i = 0; i < nDim; i++ )
    {
      star_formation_volume_min[i] = refinement_volume_min[i];
      star_formation_volume_max[i] = refinement_volume_max[i];
    }
#endif

  build_cell_buffer();
  repair_neighbors();

  /* do initial refinement */
  level = min_level;
  total_cells_per_level[min_level] = num_root_cells;
  while ( level < max_level && total_cells_per_level[level] > 0 )
    {
      cart_debug("assigning density to level %u", level );
      assign_density(level);
      cart_debug("refining level %u, num_cells_per_level = %d", level, num_cells_per_level[level] );
      modify( level, 0 );
      cart_debug("done refining level %u, created %u new cells", level, num_cells_per_level[level+1] );
      MPI_Allreduce( &num_cells_per_level[level+1], &total_cells_per_level[level+1], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
      level++;
	
      if ( local_proc_id == MASTER_NODE )
	{
	  cart_debug("level %u: %u cells", level, total_cells_per_level[level] );
	}
      
      load_balance();
    }

  if ( !buffer_enabled )
    {
      cart_debug("building cell buffer");
      build_cell_buffer();
      repair_neighbors();
    }

  for(j=min_level+1; j<=max_level; j++)
    {
      tl[j] = tl[min_level];
      aexp[j] = aexp[min_level];
    }

  run_output();

  /*
  //  Debugging parameters
  */
#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = 0;
  rt_debug.Stop = 1;
  rt_debug.Pos[0] = 29.619681;
  rt_debug.Pos[1] = 14.352527;
  rt_debug.Pos[2] = 32.880561;
#endif
#endif
}

