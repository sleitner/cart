#include "config.h"
#ifdef STAR_FORMATION

#include "auxiliary.h"
#include "control_parameter.h"
#include "imf.h"
#include "starformation_recipe.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"


#define MAX_FEEDBACK_MODELS       10
const sf_feedback_t* feedback_models[MAX_FEEDBACK_MODELS];
int num_feedback_models = 0;


double feedback_temperature_ceiling = 1.0e8;  /* Used to be called T_max_feedback; also, was a define in HART */

#ifdef BLASTWAVE_FEEDBACK
double blastwave_time = { 50.0e6 };
#endif /* BLASTWAVE_FEEDBACK */


void config_init_star_formation_feedback()
{
  int i;

  /*
  //  This must be the first call in this function
  */
  if(sf_recipe->setup_feedback != NULL) sf_recipe->setup_feedback();

  /*
  //  IMF
  */
  config_init_imf();

  /*
  // call config_init for all currently used feedback models
  */
  for(i=0; i<num_feedback_models; i++)
    {
      if(feedback_models[i]->config_init != NULL) feedback_models[i]->config_init();
    }

  /*
  //  other
  */
  control_parameter_add3(control_parameter_double,&feedback_temperature_ceiling,"fb:temperature-ceiling","feedback_temperature_ceiling","T_max_feedback","maximum gas temperature for the feedback to operate. No feedback is allowed in the gas with the temperature above this limit.");

#ifdef BLASTWAVE_FEEDBACK 
  control_parameter_add3(control_parameter_time,&blastwave_time,"blastwave-time","bw:blast-time","bw.blast_time","time before cells can cool in blastwave feedback subgrid model.");
#endif /* BLASTWAVE_FEEDBACK */
}


void config_verify_star_formation_feedback()
{
  int i;

  /*
  //  IMF
  */
  config_verify_imf();

  /*
  // call config_verify for all currently used feedback models
  */
  for(i=0; i<num_feedback_models; i++)
    {
      if(feedback_models[i]->config_verify != NULL) feedback_models[i]->config_verify();
    }

  /*
  //  other
  */
  cart_assert(feedback_temperature_ceiling > 1.0e6);

#ifdef BLASTWAVE_FEEDBACK 
  cart_assert(!(blastwave_time < 0.0));
#endif /* BLASTWAVE_FEEDBACK */
}


void add_feedback(const sf_feedback_t *ptr)
{
  int i;

  /*
  //  Do not duplicate
  */
  for(i=0; i<num_feedback_models; i++)
    {
      if(ptr == feedback_models[i]) return;
    }

  /*
  //  Is array big enough?
  */
  if(num_feedback_models == MAX_FEEDBACK_MODELS)
    {
      cart_error("Too many feedback models are implemented. Increase MAX_FEEDBACK_MODELS in starformation_feedback.c");
    }

  feedback_models[num_feedback_models++] = ptr;
}


void init_star_formation_feedback()
{
  int i;

  /*
  // call init_feedback for all currently used feedback models
  */
  for(i=0; i<num_feedback_models; i++)
    {
      if(feedback_models[i]->init_feedback != NULL) feedback_models[i]->init_feedback();
    }
}


#if defined(HYDRO) && defined(PARTICLES)

#ifdef BLASTWAVE_FEEDBACK
void init_blastwave(int icell)
{
  cell_blastwave_time(icell) =  cell_gas_density(icell) * blastwave_time;
}
#endif /* BLASTWAVE_FEEDBACK */


double dUfact;  /* must be here to simplify OpenMP directives */


void stellar_feedback(int level, int cell, int ipart, double t_next )
{
  int i;

  /*
  // call feedback for all currently used feedback models
  */
  for(i=0; i<num_feedback_models; i++)
    {
      cart_assert(feedback_models[i]->feedback != NULL);
      feedback_models[i]->feedback(level,cell,ipart,t_next);
    }
}


void setup_star_formation_feedback(int level)
{
  int i;

  dUfact = feedback_temperature_ceiling/(units->temperature*constants->wmu*(constants->gamma-1));

  /*
  // call level_setup for all currently used feedback models
  */
  for(i=0; i<num_feedback_models; i++)
    {
      if(feedback_models[i]->setup_level != NULL) feedback_models[i]->setup_level(level);
    }
}

/*
//  Currently implemented models
*/
extern sf_feedback_t sf_feedback_snII;
extern sf_feedback_t sf_feedback_snIa;
extern sf_feedback_t sf_feedback_ml;

const struct StarFormationFeedbackModels sf_feedback = 
{
  &sf_feedback_snII,
  &sf_feedback_snIa,
  &sf_feedback_ml
};

#else  /* HYDRO && PARTICLES */

const struct StarFormationFeedbackModels sf_feedback = { NULL, NULL, NULL };

#endif /* HYDRO && PARTICLES */
#endif /* STAR_FORMATION */
