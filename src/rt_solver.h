#ifndef __RT_SOLVER_H__
#define __RT_SOLVER_H__


/*
//  Usefull, RT-independent calls for diagnistics
*/
void rtSetTemUnits();
float rtTemInK(int cell);
void rtGetSobolevFactors(int cell, int level, float *len, float *vel);


#ifdef RADIATIVE_TRANSFER

#include "rt_tree.h"


#ifdef RT_VAR_SOURCE
void rtInitSource(int level);
float rtSource(int ipart);
#endif


void rtApplyCooling(int level, int num_level_cells, int *level_cells);

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
void rtGetPhotoRates(int cell, float rate[]);
void rtGetRadiationField(int cell, int n, int lxi[], float ngxi[]);
void rtGetRadiationBackground(int *nPtr, float **wlenPtr, float **ngxiPtr);

#endif  /* RADIATIVE_TRANSFER */
#endif  /* __RT_SOLVER_H__ */
