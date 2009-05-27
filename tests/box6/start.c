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

#ifdef RADIATIVE_TRANSFER
  rtSetTemUnits();
#endif

  bb[0] = 0.0;
  bb[1] = num_grid;
  bb[2] = 0.0;
  bb[3] = num_grid;
  bb[4] = 0.0;
  bb[5] = num_grid;

  sprintf(filename,"%s/out-box.%04d.bin",output_directory,(int)(auni[0]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

  extFindMaxVar(HVAR_GAS_DENSITY,&dmax,pos);

  dbb = nbin1*pow(0.5,(double)max_level);

  bb[0] = pos[0] - 0.5*dbb;
  bb[1] = pos[0] + 0.5*dbb;
  bb[2] = pos[1] - 0.5*dbb;
  bb[3] = pos[1] + 0.5*dbb;
  bb[4] = pos[2] - 0.5*dbb;
  bb[5] = pos[2] + 0.5*dbb;

  sprintf(filename,"%s/out-zoom.%04d.bin",output_directory,(int)(auni[0]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

#ifdef DEBUG_MEMORY_USE
  dmuPrintRegistryContents();
#endif

}

void init_run()
{
  int i, j;
  char filename[256], filename2[256];

#ifdef PARTICLES
  sprintf( filename, "../IC/PMcrd.DAT" );
  sprintf( filename2, "../IC/PMcrs0.DAT" );

  restart_load_balance( NULL, filename, filename2 );

  read_particles( filename, filename2, NULL, NULL, 0, NULL );
  cart_debug("read in particles");
#endif

#ifdef HYDRO
  sprintf( filename, "../IC/tr_ic.dat" );
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
  xInit[5] = (auni[min_level] < 0.0125) ? 7.0e-7/xH : 2.0e-6/xH;

  for(i=0; i<num_cells_per_level[min_level]; i++)
    {
      for(j=0; j<6; j++) cell_vars[i][RT_HVAR_OFFSET+j] = xInit[j]*cell_gas_density(i);
    }
#endif

  hydro_magic( min_level );
  hydro_eos( min_level );
#endif /* HYDRO */

  cart_debug("tl[min_level] = %f", tl[min_level] );
  cart_debug(" a[min_level] = %f", auni[min_level] );

  dtl[min_level] = 0.0;
  choose_timestep( &dtl[min_level] );

#ifdef PARTICLES
  build_mesh();
#endif /* PARTICLES */

#ifdef STARFORM
  for ( i = 0; i < nDim; i++ )
    {
      star_formation_volume_min[i] = refinement_volume_min[i];
      star_formation_volume_max[i] = refinement_volume_max[i];
    }
#endif

  if ( !buffer_enabled )
    {
      cart_debug("building cell buffer");
      build_cell_buffer();
      repair_neighbors();
    }

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

