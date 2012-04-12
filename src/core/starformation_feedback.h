#ifndef __STAR_FORMATION_FEEDBACK_H__
#define __STAR_FORMATION_FEEDBACK_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

/*
//  ATTENTION DEVELOPERS:
//  ONLY add new members at the end of the structure!!!
*/ 
struct StellarFeedback
{
  const char *name;
  void (*thermal_feedback)(int level, int cell, int ipart, double t_next);
  float (*ionizing_luminosity)(int ipart);
  float (*extra_pressure)(int cell);  /* can be NULL */
  void (*config_init)();              /* can be NULL */
  void (*config_verify)();            /* can be NULL */
  void (*init)();                     /* can be NULL */
  void (*setup)(int level);           /* can be NULL */
};


extern const struct StellarFeedback *sf_feedback;


void config_init_star_formation_feedback();
void config_verify_star_formation_feedback();

void init_star_formation_feedback();

void stellar_feedback(int level, int iter_cell, int ipart, double t_next );

void setup_star_formation_feedback(int level);

#ifdef BLASTWAVE_FEEDBACK
void init_blastwave(int icell);
#endif /* BLASTWAVE_FEEDBACK */

#endif /* STAR_FORMATION */

#endif
