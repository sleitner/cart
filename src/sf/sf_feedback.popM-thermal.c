#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "starformation_feedback.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.rad.h"
#include "models/feedback.ml.h"


void sfb_config_init()
{
  snII_config_init();
  snIa_config_init();
  ml_snl2012_config_init();
}


void sfb_config_verify()
{
  snII_config_verify();
  snIa_config_verify();
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


struct StellarFeedbackParticle sf_feedback_particle_internal = 
  {
    "popM-thermal",
    sfb_hydro_feedback,
    rad_luminosity_popM,
    NULL,
    sfb_config_init,
    sfb_config_verify,
    sfb_init,
    sfb_setup,
    NULL
  };

struct StellarFeedbackCell sf_feedback_cell_internal =
{
	NULL, NULL
};


#endif /* STAR_FORMATION */
