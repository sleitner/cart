
#ifndef __FRT_PARAMETERS_H__
#define __FRT_PARAMETERS_H__


#define IVAR_Ein		 1
#define IVAR_XHI		 2
#define IVAR_XHII		 3
#define IVAR_XHeI		 4
#define IVAR_XHeII		 5
#define IVAR_XHeIII		 6
#define IVAR_XH2		 7

#ifdef RT_8SPECIES
#define IVAR_DIM		 9
#define IVAR_XH2p		 8
#define IVAR_XHm		 9
#else
#define IVAR_DIM		 7
#endif


#define IPAR_RHO 		 1
#define IPAR_ZSOL		 2
#define IPAR_DTSH		 3
#define IPAR_DTAD		 4
#define IPAR_BIAS		 5
#define IPAR_DFAC		 6
#define IPAR_SOBL		 7
#define IPAR_NUMF		 8
#define IPAR_CELL		 9
#define IPAR_VOL		10
#define IPAR_DEB		19
#define IPAR_DIM		19


#define frtRATE_PhG2             0
#define frtRATE_PiG2             1
#define frtRATE_PhG1             2
#define frtRATE_PiG1             3
#define frtRATE_PhH1             4
#define frtRATE_PiH1             5
#define frtRATE_Ci27             6
#define frtRATE_Ci28             7
#define frtRATE_Ci29             8
#define frtRATE_Ci30             9
#define frtRATE_Ci31            10
#define frtRATE_Ci32            11
#define frtRATE_CiLW            12
#define IRATE_DIM		13

#endif
