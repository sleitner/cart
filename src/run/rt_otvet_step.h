#ifndef __RT_OTVET_STEP_H__
#define __RT_OTVET_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif

#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)


struct rtGlobalValue;

void rtStepBeginTransferOtvet(struct rtGlobalValue *abcMax);
void rtLevelUpdateTransferOtvet(int level);
void rtOtvetSplitUpdate(int level, int num_level_cells, int *level_cells);

#ifdef RT_OTVET_SAVE_FLUX
/*
//  AFTER the call to step::rtLevelUpdate the array <rt_flux> contains 
//  \vec{F}\dot\vec{n} at 6 faces of each cell for the RF component 
//  specified by <rt_flux_field> 
//  (i.e. rt_flux[...][0] = -F^x_{i-1/2}, rt_flux[...][1] = F^x_{i+1/2}, etc).
//  IMPORTANT: <rt_flux> is only defined on local cells, so you cannot 
//  differentiate it.
*/
extern int rt_flux_field;
extern float rt_flux[num_cells][num_neighbors];
#endif /* RT_OTVET_SAVE_FLUX */


#endif /* defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_OTVET_STEP_H__ */
