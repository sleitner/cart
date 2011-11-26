#ifndef __RT_TRANSFER_H__
#define __RT_TRANSFER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#ifdef RT_TRANSFER


#include <mpi.h>
#include "rt_global.h"


void rtInitRunTransfer();

/*
//  Initialization helper for the analysis mode
*/
void rtInitStepTransfer();

void rtAfterAssignDensityTransfer(int level, int num_level_cells, int *level_cells);

void rtGlobalUpdateTransfer(int top_level, MPI_Comm level_com);
void rtComputeAbsLevel(int level, int num_level_cells, int *level_cells, int freq, float **abc);


#ifdef RT_SINGLE_SOURCE
extern int rtSingleSourceLevel;
extern float rtSingleSourceValue;
extern double rtSingleSourcePos[nDim];
#endif


struct rtGlobalValue;
extern struct rtGlobalValue rtAvgRF[];
extern float rtGlobalAC[];


#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_TRANSFER_H__ */
