#ifndef __FEEDBACK_ML_H__
#define __FEEDBACK_ML_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

void ml_config_init();
void ml_config_verify();
void ml_init();
void ml_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
void ml_feedback(int level, int cell, int ipart, double t_next );
#endif /* HYDRO && PARTICLES */

/*
//  Corrections introduced by Sam Leitner in Apr 2012
*/
void ml_snl2012_config_init();
void ml_snl2012_config_verify();

#endif /* STAR_FORMATION */
#endif
