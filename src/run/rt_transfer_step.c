#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER)

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "hydro.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "rt_global.h"
#include "rt_transfer.h"
#include "sfc.h"
#include "starformation.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "../frt/frt_c.h"


extern struct rtGlobalValue rtMaxAC[];
extern struct rtGlobalValue rtAvgAC[];
extern float rtGlobalAC[];


void rtTransferSplitUpdate(int level);

#ifndef RT_VAR_SOURCE
void rtTransferUpdateUniformSource(struct rtGlobalValue *avg, MPI_Comm com);
#endif


#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
#include "rt_otvet_step.h"
#endif


int rtIsThereWork()
{
#ifdef RT_TEST
  return 1;
#endif

#ifdef RT_SINGLE_SOURCE
  return (rtSingleSourceValue > 0.0);
#endif

#if defined(PARTICLES) && defined(STARFORM)
  return 1;
#else
  return 0;
#endif /* PARTICLES && STARFORM */
}


void rtStepBeginTransfer()
{
#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
  rtStepBeginTransferOtvet(rtMaxAC);
#endif
}


void rtStepEndTransfer()
{
  int i;
  frt_real abc[rt_num_freqs];

  for(i=0; i<rt_num_freqs; i++) abc[i] = rtGlobalAC[i];

  start_time(WORK_TIMER);
  frtCall(stependtransfer)(abc);
  end_time(WORK_TIMER);
}


void rtLevelUpdateTransfer(int level)
{
  if(!rtIsThereWork()) return;

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
  rtLevelUpdateTransferOtvet(level);
#endif

  rtTransferSplitUpdate(level);
}


void rtTransferSplitUpdate(int level)
{
  int i, j, k;
  int cell;
  int num_level_cells;
  int *level_cells;
  int children[num_children];
  double new_var;
  const double factor = ((double)(1.0/(1<<nDim)));

  start_time(WORK_TIMER);

  if(level < max_level)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
#pragma omp parallel for default(none), private(i,cell,j,k,children,new_var), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i];
	  if(cell_is_refined(cell))
	    {
	      /*
	      // Average over children
	      */
	      cell_all_children(cell,children);
	      for(j=0; j<rt_num_fields; j++)
		{
		  new_var = 0.0;
		  for(k=0; k<num_children; k++)
		    {
		      new_var += cell_var(children[k],rt_field_offset+j);
		    }
		  cell_var(cell,rt_field_offset+j) = new_var*factor; 
		}
	    }
	}

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
      rtOtvetSplitUpdate(level,num_level_cells,level_cells);
#endif

      cart_free(level_cells);
    }

  end_time(WORK_TIMER);

}

#endif /* RADIATIVE_TRANSFER && RT_TRANSFER */

