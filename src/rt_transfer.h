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


void rtInitRunTransfer();
void rtStepBeginTransfer();
void rtStepEndTransfer();
void rtAfterAssignDensityTransfer(int level, int num_level_cells, int *level_cells);
void rtLevelUpdateTransfer(int level);
void rtGlobalUpdateTransfer(int top_level, MPI_Comm level_com);

void rtComputeAbsLevel(int level, int ncells, int *cells, int ifreq, float **abc);

#ifdef RT_SINGLE_SOURCE
extern int rtSingleSourceLevel;
extern float rtSingleSourceValue;
extern double rtSingleSourcePos[nDim];
#endif


struct rtGlobalValue;
extern struct rtGlobalValue rtAvgRF[];


#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_TRANSFER_H__ */
