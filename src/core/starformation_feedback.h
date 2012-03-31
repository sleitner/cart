#ifndef __STAR_FORMATION_FEEDBACK_H__
#define __STAR_FORMATION_FEEDBACK_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

struct StellarFeedback
{
  const char *name;
  float (*ionizing_luminosity)(int ipart);
  void (*hydrodynamic_feedback)(int level, int cell, int ipart, double t_next);
  void (*config_init)();           /* can be NULL */
  void (*config_verify)();         /* can be NULL */
  void (*init)();                  /* can be NULL */
  void (*setup)(int level);        /* can be NULL */
};


extern const struct StellarFeedback *sf_feedback;


void config_init_star_formation_feedback();
void config_verify_star_formation_feedback();

void set_feedback_model(const struct StellarFeedback *ptr);

void init_star_formation_feedback();

void stellar_feedback(int level, int iter_cell, int ipart, double t_next );

void setup_star_formation_feedback(int level);

#ifdef BLASTWAVE_FEEDBACK
void init_blastwave(int icell);
#endif /* BLASTWAVE_FEEDBACK */

#endif /* STAR_FORMATION */

#endif
