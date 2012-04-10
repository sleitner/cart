#ifndef __FEEDBACK_SNII_H__
#define __FEEDBACK_SNII_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

void snII_config_init();
void snII_config_verify();
void snII_init();
void snII_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
void snII_hydrodynamic_feedback(int level, int cell, int ipart, double t_next );
#endif /* HYDRO && PARTICLES */

#endif /* STARFORM */
#endif
