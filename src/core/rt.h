#ifndef __RT_H__
#define __RT_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#include <mpi.h>


void rtSetupSource(int level);
extern float (*rtSource)(int ipart);

void rtConfigInit();
void rtConfigVerify();

void rtInitRun();

/*
//  This function can be called more than once per top level step
//  to update internal RT tables.
*/ 
void rtUpdateTables();

/*
//  Initialization helper for the analysis mode
*/
void rtInitStep(double dt);

/*
//  Using somewhat less descriptive function names in the expectation
//  of the interface re-write in C++, where these functions will be
//  overloaded.
*/
void rtAfterAssignDensity1(int level);
void rtAfterAssignDensity2(int level, int num_level_cells, int *level_cells);

void rtGlobalUpdate(int top_level, MPI_Comm level_com);

/*
//  Usefull wrappers; they are slow, so should only be used for output
*/
float rtTem(int cell);
float rtDmw(int cell);
float rtDmw2(int cell);  /* this version includes the floor */
float rtUmw(int cell);
void rtGetCoolingRate(int cell, float *cooling_rate, float *heating_rate);
void rtGetPhotoRates(int cell, float *rate);
void rtGetRadiationField(int cell, int n, const float *wlen, float *ngxi);

void rtModifyTimeStep(double *dt);


#ifdef RT_TRANSFER
#define DEFINE_FRT_INTEFACE(_var_,_rawrf_) \
  frt_real _var_[FRT_DIM]; \
  frt_real _rawrf_##Buffer[rt_num_fields]; \
  frt_real *_rawrf_ = _rawrf_##Buffer
#else
#define DEFINE_FRT_INTEFACE(_var_,_rawrf_) \
  frt_real _var_[FRT_DIM]; \
  frt_real *_rawrf_##Buffer = NULL; \
  frt_real *_rawrf_ = _rawrf_##Buffer
#endif

#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_H__ */
