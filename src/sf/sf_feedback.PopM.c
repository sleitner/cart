#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "starformation_feedback.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.lum.h"
#include "models/feedback.ml.h"


void PopM_config_init()
{
  /*
  //  Type II supernova feedback
  */
  snII_config_init();

  /*
  //  Type Ia supernova feedback
  */
  snIa_config_init();

  /*
  //  ionizing luminosity
  */
  lum_config_init();

  /*
  //  mass loss
  */
  ml_config_init();
}


void PopM_config_verify()
{
  /*
  //  type II supernova feedback
  */
  snII_config_verify();

  /*
  //  type Ia supernova feedback
  */
  snIa_config_verify();

  /*
  //  ionizing luminosity
  */
  lum_config_verify();

  /*
  //  mass loss
  */
  ml_config_verify();
}


void PopM_init()
{
  snII_init();
  snIa_init();
  ml_init();
}


void PopM_setup(int level)
{
  snII_setup(level);
  snIa_setup(level);
  lum_setup(level);
  ml_setup(level);
}


#if defined(HYDRO) && defined(PARTICLES)
void PopM_hydrodynamic_feedback(int level, int cell, int ipart, double t_next )
{
  snII_hydrodynamic_feedback(level,cell,ipart,t_next);
  snIa_hydrodynamic_feedback(level,cell,ipart,t_next);
  ml_hydrodynamic_feedback(level,cell,ipart,t_next);
}
#endif /* HYDRO && PARTICLES */


struct StellarFeedback sf_feedback_internal = 
  {
    "PopM",
    lum_ionizing_luminosity,
    PopM_hydrodynamic_feedback,
    PopM_config_init,
    PopM_config_verify,
    PopM_init,
    PopM_setup
  };

#endif /* STARFORM */
