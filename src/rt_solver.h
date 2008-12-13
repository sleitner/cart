#ifndef __RT_SOLVER_H__
#define __RT_SOLVER_H__


#include "rt_tree.h"


extern float rt_XH;
extern float rt_XHe;


#ifdef RT_VAR_SOURCE
float rtSource(int ipart);
#endif


void rtApplyCooling(int level, int num_level_cells, int *level_cells);

void rtInitRun();
void rtStepBegin();
void rtStepEnd();
void rtLevelUpdate(int level, MPI_Comm local_comm);

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
float rtTemInK(int cell);
void rtGetPhotoRates(int cell, float ionRates[3], float heatRates[3]);

#endif  /* __RT_SOLVER_H__ */
