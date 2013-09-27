#ifndef __STAR_FORMATION_FEEDBACK_H__
#define __STAR_FORMATION_FEEDBACK_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

DECLARE_LEVEL_ARRAY(int,star_feedback_frequency);

/*
//  ATTENTION DEVELOPERS:
//  ONLY add new members at the end of the structure!!!
//  ONLY add new members if they are inserted in a new place in the code!!!
*/ 
struct StellarFeedbackParticle
{
  const char *name;
  void (*hydro_feedback)(int level, int cell, int ipart, double t_next);
  float (*rt_source)(int ipart);
  float (*extra_pressure)(int cell);  /* can be NULL */
  void (*config_init)();              /* can be NULL */
  void (*config_verify)();            /* can be NULL */
  void (*init)();                     /* can be NULL */
  void (*setup)(int level);           /* can be NULL */
  int (*destroy_star_particle)(int level, int icell, int ipart); /* can be NULL */
};

extern const struct StellarFeedbackParticle *sf_feedback_particle;

struct StellarFeedbackCell
{
    void (*hydro_feedback_cell)(int level, int cell, double t_next, double dt); /* can be NULL */
    void (*nonlocal_feedback_cell)(int level, int cell, double t_next, double dt); /* can be NULL */
};

extern const struct StellarFeedbackCell *sf_feedback_cell;

extern double feedback_sampling_timescale;

void config_init_star_formation_feedback();
void config_verify_star_formation_feedback();

void init_star_formation_feedback();

void setup_star_formation_feedback(int level);

#ifdef BLASTWAVE_FEEDBACK
void start_blastwave(int icell);
#endif /* BLASTWAVE_FEEDBACK */

#endif /* STAR_FORMATION */

#endif
