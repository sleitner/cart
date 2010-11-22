#ifndef __RT_SOLVER_H__
#define __RT_SOLVER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


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


#include <mpi.h>


void rtInitSource(int level);
float rtSource(int ipart);

void rtApplyCooling(int level, int num_level_cells, int *level_cells);

void rtConfigInit();
void rtConfigVerify();

void rtInitRun();
void rtStepBegin();
void rtStepEnd();
void rtLevelUpdate(int level);
void rtGlobalUpdate(int top_level, MPI_Comm level_com);

/*
//  This function can be called more than once per top level step
//  to update internal RT tables. It calls rtGlobalUpdate(...) internally,
//  so if it is used, then there is no need to call rtGlobalUpdate
*/ 
void rtUpdateTables(int top_level, MPI_Comm level_com);

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
void rtGetCoolingRate(int cell, float *cooling_rate, float *heating_rate);
void rtGetPhotoRates(int cell, float rate[]);
void rtGetRadiationField(int cell, int n, const int idxi[], float ngxi[]);
void rtGetBinIds(int n, const float wlen[], int idxi[]);
void rtGetBinWavelengths(int n, const int idxi[], float wlen[]);

void rtModifyTimeStep(double *dt);

#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_SOLVER_H__ */
