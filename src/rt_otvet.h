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


#define rtStencilSize     (2*nDim*nDim)
extern double rtStencilDist2[rtStencilSize];
extern double rtStencilDelPos[rtStencilSize][nDim];
extern double rtStencilTensor[rtStencilSize][nDim*(nDim+1)/2];


void rtInitRunTransferOtvet();

/*
//  Initialization helper for the analysis mode
*/
void rtInitStepTransferOtvet(struct rtGlobalValue *maxAC);

void rtAfterAssignDensityTransferOtvet(int level, int num_level_cells, int *level_cells);

void rtGetStencil(int level, int cell, int nb[]);

#endif /* defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_OTVET_H__ */
