#include "config.h"
#include "parallel.h"
#ifdef RADIATIVE_TRANSFER


#include "auxiliary.h"
#include "rt_global.h"
#include "timing.h"
//#include "tree.h"


void rtGlobalValueInit(struct rtGlobalValue *v, float val)
{
  int level;

  cart_assert(v != NULL);

  v->Value = val;

  for(level=min_level; level<=max_level; level++)
    {
      v->buffer[level-min_level] = val;
    }
}


void rtGlobalValueUpdate(struct rtGlobalValue *v, int level, float val)
{
  cart_assert(min_level<=level && level<=max_level);
  cart_assert(v != NULL);

  v->buffer[level-min_level] = val;
}


void rtGlobalValueCommunicate(struct rtGlobalValue *v, MPI_Op op, MPI_Comm level_com)
{
  int level;
  float val;

  cart_assert(v != NULL);

  /*
  //  Update global average
  */
  start_time(WORK_TIMER);

  /* Cannot use switch as MPI_Op may not be a simple type */
  if(op == MPI_SUM)
    {
      val = 0.0;
      for(level=min_level; level<=max_level; level++)
	{
	  val += v->buffer[level-min_level];
	}
    }
  else if(op == MPI_MAX)
    {
      val = v->buffer[0];
      for(level=min_level+1; level<=max_level; level++)
	{
	  val = max(val,v->buffer[level-min_level]);
	}
    }
  else
    {
      cart_error("Bug in rt_global.c");
    }

  end_time(WORK_TIMER);

  /*
  //  Reduce local values
  */
  start_time( COMMUNICATION_TIMER );
  MPI_Allreduce(&val,&v->Value,1,MPI_FLOAT,op,level_com);
  end_time( COMMUNICATION_TIMER );
}

#endif /* RADIATIVE_TRANSFER */

