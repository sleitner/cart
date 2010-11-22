#ifndef __RT_OTVET_H__
#define __RT_OTVET_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

struct rtGlobalValue;

void rtInitRunTransferOtvet();
void rtStepBeginTransferOtvet(struct rtGlobalValue *abcMax);
void rtAfterAssignDensityTransferOtvet(int level, int num_level_cells, int *level_cells);
void rtLevelUpdateTransferOtvet(int level);

#endif /* defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_OTVET_H__ */
