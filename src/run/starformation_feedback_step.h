#ifndef __STARFORMATION_FEEDBACK_STEP_H__
#define __STARFORMATION_FEEDBACK_STEP_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

#ifdef BLASTWAVE_FEEDBACK
void check_bwtime_precision(int level);
#endif /* BLASTWAVE_FEEDBACK */

void star_particle_feedback(int level);

#endif /* STARFORM */

#endif
