
#ifndef __FRT_INDEX_H__
#define __FRT_INDEX_H__


#ifndef __FRT_OFFSET
#error "Incorrect usage of include frt_index.h"
#endif


#define FRT_Ein					(0+__FRT_OFFSET)
#define FRT_XHI					(1+__FRT_OFFSET)
#define FRT_XHII				(2+__FRT_OFFSET)
#define FRT_XHeI				(3+__FRT_OFFSET)
#define FRT_XHeII				(4+__FRT_OFFSET)
#define FRT_XHeIII				(5+__FRT_OFFSET)
#define FRT_XH2					(6+__FRT_OFFSET)
#define FRT_XH2p				(7+__FRT_OFFSET)
#define FRT_XHm					(8+__FRT_OFFSET)

#define FRT_Density				(9+__FRT_OFFSET)
#define FRT_Metallicity				(10+__FRT_OFFSET)
#define FRT_SobolevLength			(11+__FRT_OFFSET)
#define FRT_NumericalDiffusionFactor		(12+__FRT_OFFSET)
#define FRT_OTRadiationFieldLocal		(13+__FRT_OFFSET)
#define FRT_OTRadiationFieldGlobal		(14+__FRT_OFFSET)
#define FRT_ResolutionElementVolume		(15+__FRT_OFFSET)
#define FRT_ResolutionElementSize		(16+__FRT_OFFSET)
#define FRT_SourceBias				(17+__FRT_OFFSET)
#define FRT_CoolingSuppressionFactor		(18+__FRT_OFFSET)
#define FRT_ExternalHeatingRate			(19+__FRT_OFFSET)
#define FRT_LTEFlag				(20+__FRT_OFFSET)
#define FRT_Debug				(21+__FRT_OFFSET)
#define FRT_DustToGasRatio			(22+__FRT_OFFSET)
#define FRT_Gamma				(23+__FRT_OFFSET)
#define FRT_DIM					(28)

#define FRT_RATE_HeatingHeII			(0+__FRT_OFFSET)
#define FRT_RATE_IonizationHeII			(1+__FRT_OFFSET)
#define FRT_RATE_HeatingHeI			(2+__FRT_OFFSET)
#define FRT_RATE_IonizationHeI			(3+__FRT_OFFSET)
#define FRT_RATE_HeatingHI			(4+__FRT_OFFSET)
#define FRT_RATE_IonizationHI			(5+__FRT_OFFSET)
#define FRT_RATE_IonizationO8			(6+__FRT_OFFSET)
#define FRT_RATE_DissociationLW			(13+__FRT_OFFSET)
#define FRT_RATE_DIM				(14)


#endif
