
#ifndef __FRT_PARAMETERS_H__
#define __FRT_PARAMETERS_H__


#define IVAR_Ein		 0
#define IVAR_XHI		 1
#define IVAR_XHII		 2
#define IVAR_XHeI		 3
#define IVAR_XHeII		 4
#define IVAR_XHeIII		 5
#define IVAR_XH2		 6

#ifdef RT_8SPECIES
#define IVAR_DIM		 9
#define IVAR_XH2p		 7
#define IVAR_XHm		 8
#else
#define IVAR_DIM		 7
#endif


#define IPAR_RHO		 0
#define IPAR_ZSOL		 1
#define IPAR_DTSH		 2
#define IPAR_DTAD		 3
#define IPAR_BIAS		 4
#define IPAR_DFAC		 5
#define IPAR_SOBL		 6
#define IPAR_NUMF		 7
#define IPAR_CELL		 8
#define IPAR_VOL		 9
#define IPAR_DEB		18
#define IPAR_DIM		19


#define IRATE_DIM		13


#endif
