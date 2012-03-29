#ifndef __STAR_FORMATION_FEEDBACK_H__
#define __STAR_FORMATION_FEEDBACK_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

struct StarFormationFeedback
{
  const char *name;
  void (*feedback)(int level, int cell, int ipart, double t_next);
  void (*config_init)();           /* can be NULL */
  void (*config_verify)();         /* can be NULL */
  void (*init_feedback)();         /* can be NULL */
  void (*setup_level)(int level);  /* can be NULL */
};

typedef struct StarFormationFeedback sf_feedback_t;

/*
//  This is just a convenient placeholder
*/
struct StarFormationFeedbackModels
{
  const sf_feedback_t *snII;
  const sf_feedback_t *snIa;
  const sf_feedback_t *ml;
};

extern const struct StarFormationFeedbackModels sf_feedback;

void config_init_star_formation_feedback();
void config_verify_star_formation_feedback();

void add_feedback(const sf_feedback_t *ptr);

void init_star_formation_feedback();

void stellar_feedback(int level, int iter_cell, int ipart, double t_next );

void setup_star_formation_feedback(int level);

#ifdef BLASTWAVE_FEEDBACK
void init_blastwave(int icell);
#endif /* BLASTWAVE_FEEDBACK */

#endif /* STAR_FORMATION */

#endif
