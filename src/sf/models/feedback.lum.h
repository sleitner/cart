#ifndef __FEEDBACK_LUM_H__
#define __FEEDBACK_LUM_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

void lum_config_init();
void lum_config_verify();
void lum_setup(int level);

float lum_ionizing_luminosity(int ipart);
float lum2012_ionizing_luminosity(int ipart);

#endif /* STARFORM */
#endif
