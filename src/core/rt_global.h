#ifndef __RT_GLOBAL_H__
#define __RT_GLOBAL_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#include <mpi.h>


struct rtGlobalValue
{
  float Value;
  float buffer[num_refinement_levels+1];
};


void rtGlobalValueInit(struct rtGlobalValue *v, float val); 
void rtGlobalValueUpdate(struct rtGlobalValue *v, int level, float val);
void rtGlobalValueCommunicate(struct rtGlobalValue *v, MPI_Op op, MPI_Comm level_com);


#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_GLOBAL_H__ */
