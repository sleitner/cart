#include "defs.h"

#include <math.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "hydro.h"
#include "io.h"
#include "parallel.h"
#include "particle.h"
#include "refinement_indicators.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#ifdef RT_DEBUG
#include "rt_debug.h"
#endif
#endif

#include "extra/ifrit.h"

#include "helpers.h"


extern int num_options;
extern char **options;


void run_output()
{
#ifdef HYDRO

  const int nbin1 = 256;
#ifdef RADIATIVE_TRANSFER
  int varid[] = { EXT_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, EXT_GAS_TEMPERATURE, EXT_FRACTION+RT_HVAR_OFFSET+5, EXT_CELL_LEVEL, EXT_LOCAL_PROC, HVAR_GAS_ENERGY, HVAR_MOMENTUM };
#else
  int varid[] = { EXT_FRACTION+HVAR_PRESSURE, HVAR_GAS_DENSITY, EXT_CELL_LEVEL };
#endif
  int nbin[] = { nbin1, nbin1, nbin1 };
  int zoom_level, nvars = sizeof(varid)/sizeof(int);
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

  sprintf(filename,"%s/out-box.%04d.bin",output_directory,(int)(auni[min_level]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

  extFindMaxVar(HVAR_GAS_DENSITY,&dmax,pos);

  zoom_level = max_level_now_global(MPI_COMM_WORLD);
  dbb = nbin1*pow(0.5,(double)zoom_level);

  bb[0] = pos[0] - 0.5*dbb;
  bb[1] = pos[0] + 0.5*dbb;
  bb[2] = pos[1] - 0.5*dbb;
  bb[3] = pos[1] + 0.5*dbb;
  bb[4] = pos[2] - 0.5*dbb;
  bb[5] = pos[2] + 0.5*dbb;

  sprintf(filename,"%s/out-zoom.%04d.bin",output_directory,(int)(auni[min_level]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

#endif  /* HYDRO */
}


void init_run()
{
  int i, level;
  const char *rootname;
  char type[2];
  int dc_off = 0;

  /*
  //  Where do we get the root name? Use options for now
  */
  if(local_proc_id == MASTER_NODE)
    {
      if(num_options < 1)
	{
          cart_error("An option -root=<name> is required, where <name> is the root name for a set of GIC input files.\n");
        }
    }

  for(i=0; i<num_options; i++)
    {
      /*
      //  Root name for data files
      */
      rootname = check_option1(options[i],"root",NULL);
      if(rootname != NULL) continue;

      /*
      //  Switch off the DC mode
      */
      if(check_option0(options[i],"no-dc") != NULL)
	{
	  dc_off = 1;
	  continue;
	}

      cart_error("Unrecognized option: %s",options[i]);
    }

  type[1] = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  /*
  //  Begin with load balancing 
  */
#ifdef HYDRO
  type[0] = 'D';
#else
  type[0] = 'M';
#endif

  gicBalanceLoad(rootname,type);
  gicReadParticleData(rootname,type);

  cart_debug("read in particles");

#ifdef HYDRO

  type[0] = 'B';
  gicReadGasData(rootname,type);
  cart_debug("read in gas");

  hydro_magic(min_level);
  hydro_eos(min_level);

#endif /* HYDRO */

  if(dc_off) cosmology_set(DeltaDC,0.0);

  cart_debug("tl[min_level] = %f", tl[min_level] );
  cart_debug("au[min_level] = %f", auni[min_level] );
  cart_debug("ab[min_level] = %f", abox[min_level] );
  cart_debug("DC mode = %f", cosmology->DeltaDC );

  for(level=min_level+1; level<=max_level; level++)
    {
      tl[level] = tl[min_level];
      auni[level] = auni[min_level];
      abox[level] = abox[min_level];
    }

  dtl[min_level] = 0.0;
  choose_timestep( &dtl[min_level] );

  for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
    {
      particle_t[i] = tl[min_level];
      particle_dt[i] = dtl[min_level];
    }

#ifdef STARFORM
  for(i=0; i<nDim; i++)
    {
      star_formation_volume_min[i] = refinement_volume_min[i];
      star_formation_volume_max[i] = refinement_volume_max[i];
    }
#endif

  /*
  //  Debugging parameters
  */
#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = 0;
  rt_debug.Stop = 0;
  rt_debug.Pos[0] = 12.636;
  rt_debug.Pos[1] = 14.356;
  rt_debug.Pos[2] = 14.56;
#endif
#endif
}

