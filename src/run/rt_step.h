#ifndef __RT_STEP_H__
#define __RT_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


void rtApplyCooling(int level, int num_level_cells, int *level_cells);

void rtStepBegin();
void rtStepEnd();
void rtLevelUpdate(int level);

void rtModifyTimeStep(double *dt);

#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_STEP_H__ */
