#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "tree.h"
#include "units.h"
#include "starformation_feedback.h"
#include "models/onestarfits.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.ml.h"
#include "models/feedback.winds.h"
#include "models/feedback.rapSR.h"

#include "models/feedback.irtrapping.h"
#include "models/feedback.kinetic.h"
#include "models/feedback.rad.h"

#include "times.h"
#include "starformation.h"
#include "particle.h"

extern double rapSR_boost;
extern double wind_momentum_boost;
extern double tauIR_boost;

void sfb_config_init()
{
  snII_config_init();
  snIa_config_init();
  ml_snl2012_config_init();
  wind_config_init();
  rapSR_config_init();

  trapIR_config_init();
  kfb_config_init();
}


void sfb_config_verify()
{
  snII_config_verify();
  snIa_config_verify();
  ml_snl2012_config_verify();
  wind_config_verify();
  rapSR_config_verify();

  trapIR_config_verify();
  kfb_config_verify();
}


void sfb_init()
{
  snII_init();
  snIa_init();
  ml_init();
  wind_init();
  rapSR_init(); 

/*   trapIR_init();  */
  kfb_init();
}

void sfb_setup(int level)
{
    snII_setup(level);
    snIa_setup(level);
    rad_setup(level);
    ml_setup(level);
    wind_setup(level);
    rapSR_setup(level);

    trapIR_setup(level);
/*     kfb_setup(level); */
}


#if defined(HYDRO) && defined(PARTICLES)
void sfb_hydro_feedback(int level, int cell, int ipart, double t_next )
{

    snII_thermal_feedback(level,cell,ipart,t_next);
    snII_kinetic_feedback(level,cell,ipart,t_next); 
    ml_feedback(level,cell,ipart,t_next);

    if(rapSR_boost > 0){
	rapSR_kick(level,cell,ipart,t_next);    
    }

    if(wind_momentum_boost > 0){
	stellar_wind_kick(level,cell,ipart,t_next);
    }

    snIa_thermal_feedback(level,cell,ipart,t_next); 

}

extern double sf_min_gas_number_density;
void sfb_hydro_feedback_cell(int level, int cell, double t_next, double dt )
{
    if( cell_gas_density(cell)*units->number_density*constants->XH 
	> sf_min_gas_number_density){ /* note this parameter is active for any SFP */
	if(tauIR_boost>0){
	    cell_trapIR(level, cell, t_next, dt); 
	} 
    }
}
#endif /* HYDRO && PARTICLES */

struct StellarFeedbackParticle sf_feedback_particle_internal =
  {
    "popM-cluster",
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
    sfb_hydro_feedback_cell 
};

#endif /* STARFORM */

