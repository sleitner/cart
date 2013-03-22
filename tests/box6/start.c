#include "config.h"

#include "../../cosmology/hart/hart_init.h"

#include "rt_debug.h"
#include "tree.h"


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
  hart_init();
}

