#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

/* NG: definitions of pc and GMsun are from http://ssd.jpl.nasa.gov/?constants */

#define	year	(365.25*86400)	/* Julian year in seconds */
#define myr	(1e6*year)
#define gyr	(1e9*year)
#define pc	(3.0856775813e18)	/* in cm */
#define kpc	(1e3*pc)
#define mpc	(1e6*pc)

#define G	(6.67259e-8)
#define clight	(2.99792458e10)	/* cm/s */

#define Msun	(1.32712440018e26/G)	/* in grams */

#define Y_p	(0.24)	/* He mass fraction */
#define wmu	(4.0/(8.0-5.0*Y_p)) /* mol weight */
#define wmu_e	(1.0/(1.0-0.5*Y_p))

#define T_CMB0	(2.726)

#define Zsolar	(0.0199)	/* solar metallicity */

#define keV	(8.6170868e-8)

#define small	(1e-30)

#endif
