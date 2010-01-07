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
void rtLevelUpdateTransfer(int level, MPI_Comm local_comm);

void rtComputeAbsLevel(int ncells, int *cells, int ifreq, float **abc);

#ifdef RT_SINGLE_SOURCE
extern float rtSingleSourceVal;
extern double rtSingleSourcePos[nDim];
#endif

#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_TRANSFER_H__ */
