#include "config.h"

#include <math.h>
#include <stdio.h>

#include "io.h"
#include "timestep.h"
#include "tree.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#ifdef RT_DEBUG
#include "rt_debug.h"
#endif
#endif


#include "gic_init.h"


void run_output()
{
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

