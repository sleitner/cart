#ifndef __FEEDBACK_ML_H__
#define __FEEDBACK_ML_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

void ml_config_init();
void ml_config_verify();
void ml_init();
void ml_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
void ml_hydrodynamic_feedback(int level, int cell, int ipart, double t_next );
#endif /* HYDRO && PARTICLES */

#endif /* STARFORM */
#endif
