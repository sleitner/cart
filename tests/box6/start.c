#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "density.h"
#include "gravity.h"
#include "hydro.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "sfc.h"
#include "starformation.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "rt_debug.h"
#include "rt.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


#include "extra/ifrit.h"


void run_output()
{
  /*
  //  Debugging parameters
  */
#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = 0;
  rt_debug.Stop = 2;
  cell_position_double(8628465,rt_debug.Pos);
#endif
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

  units_reset();
  units_update(min_level);

#ifdef RADIATIVE_TRANSFER
  float xInit[6];

#ifdef COSMOLOGY
  xInit[1] = 1.2e-5*sqrt(cosmology->OmegaM)/(cosmology->OmegaB*cosmology->h);
#else
 #error "COSMOLOGY must be specified."
#endif

  xInit[0] = constants->XH - xInit[1];
  xInit[4] = xInit[3] = 1.0e-10;
  xInit[2] = constants->XHe - xInit[3] - xInit[4];
  xInit[5] = (auni[min_level] < 0.0125) ? 7.0e-7*constants->XH : 2.0e-6*constants->XH;

  for(i=0; i<num_cells_per_level[min_level]; i++)
    {
      for(j=0; j<6; j++) cell_vars[i][RT_HVAR_OFFSET+j] = xInit[j]*cell_gas_density(i);
    }
#endif

  hydro_magic( min_level );
  hydro_eos( min_level );
#endif /* HYDRO */

  cart_debug("tl[min_level] = %f", tl[min_level] );
#ifdef COSMOLOGY
  cart_debug(" a[min_level] = %f", auni[min_level] );
#endif

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
}

