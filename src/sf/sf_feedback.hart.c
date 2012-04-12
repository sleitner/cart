#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "starformation_feedback.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.rad.h"
#include "models/feedback.ml.h"


void sfb_config_init()
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
  //  mass loss
  */
  ml_config_init();
}


void sfb_config_verify()
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
  //  mass loss
  */
  ml_config_verify();
}


void sfb_init()
{
  snII_init();
  snIa_init();
  ml_init();
}


void sfb_setup(int level)
{
  snII_setup(level);
  snIa_setup(level);
  rad_setup(level);
  ml_setup(level);
}


#if defined(HYDRO) && defined(PARTICLES)
void sfb_hydro_feedback(int level, int cell, int ipart, double t_next )
{
  snII_thermal_feedback(level,cell,ipart,t_next);
  snIa_thermal_feedback(level,cell,ipart,t_next);
  ml_feedback(level,cell,ipart,t_next);
}
#endif /* HYDRO && PARTICLES */


struct StellarFeedback sf_feedback_internal = 
  {
    "hart",
    sfb_hydro_feedback,
    rad_luminosity_hart,
    NULL,
    sfb_config_init,
    sfb_config_verify,
    sfb_init,
    sfb_setup
  };

#endif /* STARFORM */
