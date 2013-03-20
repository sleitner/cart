#include "config.h"
#ifdef STAR_FORMATION

#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "imf.h"
#include "starformation_recipe.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"

extern struct StellarFeedbackParticle sf_feedback_particle_internal;
const struct StellarFeedbackParticle *sf_feedback_particle = &sf_feedback_particle_internal;

extern struct StellarFeedbackCell sf_feedback_cell_internal;
const struct StellarFeedbackCell *sf_feedback_cell = &sf_feedback_cell_internal;

double feedback_temperature_ceiling = 1.0e8;  /* Used to be called T_max_feedback; also, was a define in HART */
double feedback_speed_time_ceiling = 1e4;  /* do not accelerate flows so that time-steps drop below~1e4yr (1e3km/s at 10pc) */

#ifdef BLASTWAVE_FEEDBACK
double blastwave_time = { 50.0e6 };
#endif /* BLASTWAVE_FEEDBACK */


void control_parameter_list_feedback(FILE *stream, const void *ptr)
{
  fprintf(stream,"<%s>",sf_feedback_particle->name);
}


void config_init_star_formation_feedback()
{
  static char *ptr;
  ControlParameterOps r = { NULL, control_parameter_list_feedback };

  /*
  // call config_init for the current feedback model
  */
  if(sf_feedback_particle->config_init != NULL) sf_feedback_particle->config_init();

  /*
  //  other
  */
  ptr = cart_alloc(char,strlen(sf_feedback_particle_internal.name)+1);
  strcpy(ptr,sf_feedback_particle_internal.name);
  control_parameter_add(r,ptr,"sf:feedback","a feedback model for star formation. This parameter is for listing only, and must be set with SF_FEEDBACK define in defs.h. See /src/sf for available feedback models.");

  control_parameter_add2(control_parameter_time,&feedback_speed_time_ceiling,"fb:time-ceiling","feedback_speed_time_ceiling","minimum cell crossing time feedback can contribute to. Kinetic feedback doesn't add to gas speed above this limit.");

  control_parameter_add3(control_parameter_double,&feedback_temperature_ceiling,"fb:temperature-ceiling","feedback_temperature_ceiling","T_max_feedback","maximum gas temperature for the feedback to operate. No feedback is allowed in the gas with the temperature above this limit.");


#ifdef BLASTWAVE_FEEDBACK 
  control_parameter_add3(control_parameter_time,&blastwave_time,"blastwave-time","bw:blast-time","bw.blast_time","time before cells can cool in blastwave feedback subgrid model.");
#endif /* BLASTWAVE_FEEDBACK */
}


#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)

void config_verify_star_formation_feedback()
{
  char feedback_internal_name[99];
#ifdef SF_FEEDBACK
  const char *feedback_external_name = to_string(SF_FEEDBACK);
#else
  const char *feedback_external_name = "";
#endif

  cart_assert(sf_feedback_particle_internal.name != NULL);
  cart_assert(sf_feedback_particle_internal.rt_source != NULL);
  cart_assert(sf_feedback_particle_internal.hydro_feedback != NULL);

  sprintf(feedback_internal_name,"<%s>",sf_feedback_particle_internal.name);
  if(strcmp("<custom>",feedback_external_name)!=0 && strcmp(feedback_internal_name,feedback_external_name)!=0)
    {
      cart_error("Misconfiguration: the internal SF feedback name (%s) does not match the name set in defs.h (%s)",feedback_internal_name,feedback_external_name);
    }

  /*
  // call config_verify for the current feedback model
  */
  if(sf_feedback_particle->config_verify != NULL) sf_feedback_particle->config_verify();

  /*
  //  other
  */
  VERIFY(fb:temperature-ceiling, feedback_temperature_ceiling > 1.0e6 );
  VERIFY(fb:time-ceiling, feedback_speed_time_ceiling < 1.0e6 && feedback_speed_time_ceiling >0);

#ifdef BLASTWAVE_FEEDBACK 
  VERIFY(blastwave-time, !(blastwave_time < 0.0) );
#endif /* BLASTWAVE_FEEDBACK */
}


void init_star_formation_feedback()
{
  /*
  // call init for the current feedback model
  */
  if(sf_feedback_particle->init != NULL) sf_feedback_particle->init();
}


#if defined(HYDRO) && defined(PARTICLES)

#ifdef BLASTWAVE_FEEDBACK
void init_blastwave(int icell)
{
  cell_blastwave_time(icell) =  cell_gas_density(icell) * blastwave_time;
}
#endif /* BLASTWAVE_FEEDBACK */


double dUfact;  /* must be here to simplify OpenMP directives */
double dvfact;  

void stellar_feedback_particle(int level, int cell, int ipart, double t_next )
{
  /*
  // call feedback for the current feedback model for each particle
  */
  sf_feedback_particle->hydro_feedback(level,cell,ipart,t_next);
}

void stellar_feedback_cell(int level, int cell, double t_next, double dt )
{
  /*
  // call feedback for the current feedback model for each cell
  */
    sf_feedback_cell->hydro_feedback_cell(level,cell,t_next,dt);
}

/* void stellar_destruction(int level, int cell,  int ipart, int *icheck ) { */
/*   /\* */
/*   // call particle destruction for the current feedback model  */
/*   *\/ */
/*     sf_feedback_particle->destroy_star_particle(level,cell,ipart,icheck); */
/* } */

void setup_star_formation_feedback(int level)
{
  dUfact = feedback_temperature_ceiling/(units->temperature*constants->wmu*(constants->gamma-1));
  /* speed ceiling -- multiply by density for momentum ceiling */
  dvfact = cell_size[level] / (feedback_speed_time_ceiling*constants->yr/units->time);

  /*
  // call setup for the current feedback model
  */
  if(sf_feedback_particle->setup != NULL) sf_feedback_particle->setup(level);
}

#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
