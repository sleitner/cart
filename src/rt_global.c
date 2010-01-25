#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER)

#include "parallel.h"
#include "rt_global.h"
#include "timing.h"
#include "tree.h"


#ifdef RT_VAR_OT_FIELD

struct rtGlobalAverageData rtGlobals[2];
int rtNumGlobals = 2;

#else

struct rtGlobalAverageData rtGlobals[1];
int rtNumGlobals = 1;

#endif


void rtGlobalAverageInit(int n, struct rtGlobalAverageData *out)
{
  int i, level;

  start_time(WORK_TIMER);

  for(i=0; i<n; i++)
    {
      out[i].Value = 0.0;
      for(level=min_level; level<=max_level; level++)
	{
	  out[i].LocalLevelSum[level-min_level] = out[i].GlobalLevelSum[level-min_level] = 0.0;
	}
    }

  end_time(WORK_TIMER);
}


void rtGlobalAverageUpdate(int level, int n, struct rtGlobalAverageData *out, MPI_Comm local_comm)
{
  int i, lev;

  for(i=0; i<n; i++)
    {
      /*
      //  Reduce local values
      */
      start_time( COMMUNICATION_TIMER );
      MPI_Allreduce(out[i].LocalLevelSum+(level-min_level),out[i].GlobalLevelSum+(level-min_level),1,MPI_FLOAT,MPI_SUM,local_comm);
      end_time( COMMUNICATION_TIMER );

      /*
      //  Update global average
      */
      start_time(WORK_TIMER);
      out[i].Value = 0.0;
      for(lev=min_level; lev<=max_level; lev++)
	{
	  out[i].Value += out[i].GlobalLevelSum[lev-min_level];
	}
      out[i].Value /= num_root_cells;
      end_time(WORK_TIMER);
    }
}

#endif /* RADIATIVE_TRANSFER && RT_TRANSFER */

