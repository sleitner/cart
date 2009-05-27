#include "defs.h"

#include <math.h>

#include "io.h"
#include "timestep.h"
#include "tree.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#ifdef RT_DEBUG
#include "rt_debug.h"
#endif
#endif

#include "extra/ifrit.h"


#include "gic_init.h"


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
  gic_init();
  
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

