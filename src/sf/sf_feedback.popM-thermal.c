#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "starformation_feedback.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.lum.h"
#include "models/feedback.ml.h"


void sfb_config_init()
{
  snII_config_init();
  snIa_config_init();
  lum_config_init();
  ml_snl2012_config_init();
}


void sfb_config_verify()
{
  snII_config_verify();
  snIa_config_verify();
  lum_config_verify();
  ml_snl2012_config_verify();
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
  lum_setup(level);
  ml_setup(level);
}


#if defined(HYDRO) && defined(PARTICLES)
void sfb_thermal_feedback(int level, int cell, int ipart, double t_next )
{
  snII_thermal_feedback(level,cell,ipart,t_next);
  snIa_thermal_feedback(level,cell,ipart,t_next);
  ml_thermal_feedback(level,cell,ipart,t_next);
}
#endif /* HYDRO && PARTICLES */


struct StellarFeedback sf_feedback_internal = 
  {
    "popM-thermal",
    sfb_thermal_feedback,
    lum_ionizing_luminosity_popM,
    NULL,
    sfb_config_init,
    sfb_config_verify,
    sfb_init,
    sfb_setup
  };

#endif /* STARFORM */
