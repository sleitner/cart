#ifndef __STARFORMATION_FEEDBACK_H__
#define __STARFORMATION_FEEDBACK_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

void config_init_star_formation_feedback();
void config_verify_star_formation_feedback();

void init_star_formation_feedback();

void stellar_feedback(int level, int iter_cell, int ipart, double t_next );

void setup_star_formation_feedback(int level);

#ifdef BLASTWAVE_FEEDBACK
void init_blastwave(int icell);
void check_bwtime_precision(int level);
#endif /* BLASTWAVE_FEEDBACK */

void nonthermal_particle_feedback(int level, int cell, int ipart, double t_next); 

#endif /* STARFORM */

#endif
