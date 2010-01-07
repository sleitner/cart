#ifndef __RT_SOLVER_H__
#define __RT_SOLVER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>

/*
//  Usefull, RT-independent calls for diagnistics
*/
#ifdef HYDRO
void rtGetSobolevFactors(int cell, int level, float *len, float *vel);
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


void rtInitSource(int level);
float rtSource(int ipart);

void rtApplyCooling(int level, int num_level_cells, int *level_cells);

void rtConfigInit();
void rtConfigVerify();

void rtInitRun();
void rtStepBegin();
void rtStepEnd();
void rtLevelUpdate(int level, MPI_Comm local_comm);

/*
//  This function can be called more than once per top level step
//  to update internal RT tables.
*/ 
void rtUpdateTables();

/*
//  Using somewhat less descriptive function names in the expectation
//  of the interface re-write in C++, where these functions will be
//  overloaded.
*/
void rtAfterAssignDensity1(int level);
void rtAfterAssignDensity2(int level, int num_level_cells, int *level_cells);

/*
//  Usefull wrappers; they are slow, so should only be used for output
*/
float rtTem(int cell);
float rtDustToGas(int cell);
void rtGetPhotoRates(int cell, float rate[]);
void rtGetRadiationField(int cell, int n, const int idxi[], float ngxi[]);
void rtGetBinIds(int n, const float wlen[], int idxi[]);
void rtGetBinWavelengths(int n, const int idxi[], float wlen[]);

#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_SOLVER_H__ */
