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

#ifdef RT_OTVET_SAVE_FLUX
extern int rt_flux_frequency;
extern float rt_flux[num_cells][nDim];
#endif /* RT_OTVET_SAVE_FLUX */


#endif /* defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_OTVET_STEP_H__ */
