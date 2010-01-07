#ifndef __POWER_H__
#define __POWER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


void compute_power_spectrum( char *filename, int power_type );

#define POWER_TYPE_TOTAL	0
#define POWER_TYPE_DARK		1
#define POWER_TYPE_BARYONS	2
#define POWER_TYPE_GAS		3
#define POWER_TYPE_STARS	4

#endif
