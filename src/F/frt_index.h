
#ifndef __FRT_INDEX_H__
#define __FRT_INDEX_H__


#ifndef frtOFFSET
#error "Incorrect usage of include frt_index.h"
#endif


#define frtVAR_Ein		 (0+frtOFFSET)
#define frtVAR_XHI		 (1+frtOFFSET)
#define frtVAR_XHII		 (2+frtOFFSET)
#define frtVAR_XHeI		 (3+frtOFFSET)
#define frtVAR_XHeII		 (4+frtOFFSET)
#define frtVAR_XHeIII		 (5+frtOFFSET)
#define frtVAR_XH2		 (6+frtOFFSET)

#ifdef RT_8SPECIES
#define frtVAR_DIM		 (9+frtOFFSET)
#define frtVAR_XH2p		 (7+frtOFFSET)
#define frtVAR_XHm		 (8+frtOFFSET)
#else
#define frtVAR_DIM		 (7+frtOFFSET)
#endif


#define frtPAR_RHO		 (0+frtOFFSET)
#define frtPAR_ZSOL		 (1+frtOFFSET)
#define frtPAR_DTSH		 (2+frtOFFSET)
#define frtPAR_DTAD		 (3+frtOFFSET)
#define frtPAR_BIAS		 (4+frtOFFSET)
#define frtPAR_DFAC		 (5+frtOFFSET)
#define frtPAR_SOBL		 (6+frtOFFSET)
#define frtPAR_NUMF		 (7+frtOFFSET)
#define frtPAR_CELL		 (8+frtOFFSET)
#define frtPAR_VOL		 (9+frtOFFSET)
#define frtPAR_DEB		(18+frtOFFSET)
#define frtPAR_DIM		(19+frtOFFSET)


#define frtRATE_PhG2             (0+frtOFFSET)
#define frtRATE_PiG2             (1+frtOFFSET)
#define frtRATE_PhG1             (2+frtOFFSET)
#define frtRATE_PiG1             (3+frtOFFSET)
#define frtRATE_PhH1             (4+frtOFFSET)
#define frtRATE_PiH1             (5+frtOFFSET)
#define frtRATE_CiLW            (12+frtOFFSET)
#define frtRATE_DIM		(13+frtOFFSET)


#endif
