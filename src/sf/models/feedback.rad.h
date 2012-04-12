#ifndef __FEEDBACK_RAD_H__
#define __FEEDBACK_RAD_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

void rad_setup(int level);

float rad_luminosity_hart(int ipart);
float rad_luminosity_popM(int ipart);

#endif /* STARFORM */
#endif
