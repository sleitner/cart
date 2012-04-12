#ifndef __FEEDBACK_LUM_H__
#define __FEEDBACK_LUM_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

void lum_config_init();
void lum_config_verify();
void lum_setup(int level);

float lum_ionizing_luminosity_hart(int ipart);
float lum_ionizing_luminosity_popM(int ipart);

#endif /* STARFORM */
#endif
