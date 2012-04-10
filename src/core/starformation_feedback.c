#include "config.h"
#ifdef STAR_FORMATION

#include "auxiliary.h"
#include "control_parameter.h"
#include "imf.h"
#include "starformation_recipe.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"


const struct StellarFeedback *sf_feedback = NULL;


double feedback_temperature_ceiling = 1.0e8;  /* Used to be called T_max_feedback; also, was a define in HART */

#ifdef BLASTWAVE_FEEDBACK
double blastwave_time = { 50.0e6 };
#endif /* BLASTWAVE_FEEDBACK */


void config_init_star_formation_feedback()
{
  cart_assert(sf_feedback != NULL);

  /*
  //  IMF
  */
  config_init_imf();

  /*
  // call config_init for the current feedback model
  */
  if(sf_feedback->config_init != NULL) sf_feedback->config_init();

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
  /*
  //  IMF
  */
  config_verify_imf();

  /*
  // call config_verify for the current feedback model
  */
  if(sf_feedback->config_verify != NULL) sf_feedback->config_verify();

  /*
  //  other
  */
  VERIFY(fb:temperature-ceiling, feedback_temperature_ceiling > 1.0e6 );

#ifdef BLASTWAVE_FEEDBACK 
  VERIFY(blastwave-time, !(blastwave_time < 0.0) );
#endif /* BLASTWAVE_FEEDBACK */
}


void set_feedback_model(const struct StellarFeedback *ptr)
{
  /*
  //  This functions allows to replace the default feedback model
  */
  sf_feedback = ptr;

  cart_assert(sf_feedback);
  cart_assert(sf_feedback->ionizing_luminosity != NULL);
  cart_assert(sf_feedback->hydrodynamic_feedback != NULL);
}


void init_star_formation_feedback()
{
  /*
  // call init for the current feedback model
  */
  if(sf_feedback->init != NULL) sf_feedback->init();
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
  /*
  // call hydrodynamic_feedback for the current feedback model
  */
  sf_feedback->hydrodynamic_feedback(level,cell,ipart,t_next);
}


void setup_star_formation_feedback(int level)
{
  dUfact = feedback_temperature_ceiling/(units->temperature*constants->wmu*(constants->gamma-1));

  /*
  // call setup for the current feedback model
  */
  if(sf_feedback->setup != NULL) sf_feedback->setup(level);
}

#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
