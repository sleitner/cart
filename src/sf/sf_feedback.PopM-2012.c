#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "starformation_feedback.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.lum.h"
#include "models/feedback.ml2012.h"


void PopM2012_config_init()
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
  ml2012_config_init();
}


void PopM2012_config_verify()
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
  ml2012_config_verify();
}


void PopM2012_init()
{
  snII_init();
  snIa_init();
  ml2012_init();
}


void PopM2012_setup(int level)
{
  snII_setup(level);
  snIa_setup(level);
  lum_setup(level);
  ml2012_setup(level);
}


#if defined(HYDRO) && defined(PARTICLES)
void PopM2012_hydrodynamic_feedback(int level, int cell, int ipart, double t_next )
{
  snII_hydrodynamic_feedback(level,cell,ipart,t_next);
  snIa_hydrodynamic_feedback(level,cell,ipart,t_next);
  ml2012_hydrodynamic_feedback(level,cell,ipart,t_next);
}
#endif /* HYDRO && PARTICLES */


struct StellarFeedback sf_feedback_internal = 
  {
    "PopM-2012",
    lum2012_ionizing_luminosity,
    PopM2012_hydrodynamic_feedback,
    PopM2012_config_init,
    PopM2012_config_verify,
    PopM2012_init,
    PopM2012_setup
  };

#endif /* STARFORM */
