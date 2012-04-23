#ifndef __FEEDBACK_SNIa_H__
#define __FEEDBACK_SNIa_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

void snIa_config_init();
void snIa_config_verify();
void snIa_init();
void snIa_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
void snIa_thermal_feedback(int level, int cell, int ipart, double t_next );
#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
#endif
