#ifndef __STAR_FORMATION_FEEDBACK_STEP_H__
#define __STAR_FORMATION_FEEDBACK_STEP_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

#ifdef BLASTWAVE_FEEDBACK
void check_bwtime_precision(int level);
#endif /* BLASTWAVE_FEEDBACK */

void star_particle_feedback(int level);

#endif /* STAR_FORMATION */

#endif
