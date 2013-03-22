#ifndef __FEEDBACK_AGB_H__
#define __FEEDBACK_AGB_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

void AGB_config_init();
void AGB_config_verify();
void AGB_init();
void AGB_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
void AGB_feedback(int level, int cell, int ipart, double t_next );
#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
#endif
