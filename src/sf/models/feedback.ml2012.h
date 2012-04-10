#ifndef __FEEDBACK_ML2012_H__
#define __FEEDBACK_ML2012_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "feedback.ml.h"


#ifdef STARFORM

void ml2012_config_init();
void ml2012_config_verify();
void ml2012_init();
void ml2012_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
/*
//  Used with a #define for efficiency
*/
#define ml2012_hydrodynamic_feedback  ml_hydrodynamic_feedback
#endif /* HYDRO && PARTICLES */

#endif /* STARFORM */
#endif
