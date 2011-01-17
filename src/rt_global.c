#include "config.h"
#include "parallel.h"
#ifdef RADIATIVE_TRANSFER


#include "auxiliary.h"
#include "rt_global.h"
#include "timing.h"
//#include "tree.h"


void rtGlobalValueInit(struct rtGlobalValue *v, float val)
{
  int l;

  cart_assert(v != NULL);

  v->Value = val;

  for(l=0; l<=max_level-min_level; l++)
    {
      v->lb[l] = v->gb[l] = val;
    }
}


void rtGlobalValueChange(struct rtGlobalValue *v, int level, float val)
{
  cart_assert(min_level<=level && level<=max_level);
  cart_assert(v != NULL);

  v->lb[level-min_level] = val;
}


void rtGlobalValueUpdate(struct rtGlobalValue *v, int level1, int level2, MPI_Op op, MPI_Comm level_com)
{
  int l;

  cart_assert(min_level<=level1 && level1<=level2 && level2<=max_level);
  cart_assert(v != NULL);

  /*
  //  Reduce local values
  */
  start_time( COMMUNICATION_TIMER );
  MPI_Allreduce(v->lb+level1-min_level,v->gb+level1-min_level,level2-level1+1,MPI_FLOAT,op,level_com);
  end_time( COMMUNICATION_TIMER );

  /*
  //  Update global average
  */
  start_time(WORK_TIMER);

  if(op == MPI_SUM)
    {
      v->Value = 0.0;
      for(l=0; l<=max_level-min_level; l++)
	{
	  v->Value += v->gb[l];
	}
    }
  else if(op == MPI_MAX)
    {
      v->Value = v->gb[0];
      for(l=1; l<=max_level-min_level; l++)
	{
	  v->Value = max(v->Value,v->gb[l]);
	}
    }
  else
    {
      cart_error("Bug in rt_global.c");
    }

  end_time(WORK_TIMER);
}

#endif /* RADIATIVE_TRANSFER */

