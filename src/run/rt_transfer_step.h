#ifndef __RT_TRANSFER_STEP_H__
#define __RT_TRANSFER_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif

#ifdef RT_TRANSFER


void rtStepBeginTransfer();
void rtStepEndTransfer();
void rtLevelUpdateTransfer(int level);

#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_TRANSFER_STEP_H__ */
