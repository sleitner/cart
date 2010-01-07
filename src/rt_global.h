#ifndef __RT_GLOBAL_H__
#define __RT_GLOBAL_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#ifdef RT_TRANSFER

#include <mpi.h>


#define RT_SOURCE_AVG      0
#define RT_OT_FIELD_AVG    1


struct rtGlobalAverageData
{
  float Value;
  float LocalLevelSum[num_refinement_levels+1];
  float GlobalLevelSum[num_refinement_levels+1];
};

extern struct rtGlobalAverageData rtGlobals[];
extern int rtNumGlobals;


void rtGlobalAverageInit(int n, struct rtGlobalAverageData *out); 
void rtGlobalAverageUpdate(int level, int n, struct rtGlobalAverageData *out, MPI_Comm local_comm); 

#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_GLOBAL_H__ */
