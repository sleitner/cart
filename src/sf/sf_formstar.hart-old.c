#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>
#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_formstar.h"
#include "rand.h"
#include "tree.h"
#include "units.h"

#include "models/form_star.hart.h"

extern double sf_sampling_timescale; /* just check sf_timescale is >2x longer */
extern double sf_timescale;                 /* in yrs; used to be called tau_SF, did not exist in HART */
extern double sf_min_stellar_particle_mass;   /* in Msun; used to be called dm_star_min */
extern double mstar_min;
extern double dt_SF;
void star_form_config_init()
{
    formstar_hart_config_init();
}


void star_form_config_verify()
{
    formstar_hart_config_verify();
}


void star_form_setup(int level)
{
    formstar_hart_setup(level);
}

void star_form_particles(int level, int icell, double dtl, double dt_eff, float sfr){
    int i;
    double mstar;
    double P_SF, P_mass;
    
    /* probability of forming a star is Poisson with <t> = dt_SF */
    P_SF = exp( -dt_eff / dt_SF );

    mstar = sfr*dt_SF*cell_volume[level];
    
    
    if(mstar < mstar_min) P_mass = 1 - mstar/mstar_min; else P_mass = 0;
    
    /* randomly generate particle on timescale dt_SF */
    if(cart_rand() > P_SF+P_mass-P_SF*P_mass)
        mstar = MAX(mstar_min,mstar);
    else 
        mstar = 0.0;
    
    if ( mstar > 0.0 ) {
        
        /* create the new star */
        create_star_particle( icell, mstar, dtl, STAR_TYPE_NORMAL );
        
#ifdef LOG_STAR_CREATION      
        log_star_creation( icell, particle_mass[ipart], FILE_RECORD);
#endif
        
    }
}

struct FormStar sf_formstar_internal = 
  {
    "hart-old\0",
    star_form_particles,
    star_form_config_init,
    star_form_config_verify,
    NULL, 
    star_form_setup
  };

#endif /* HYDRO && STAR_FORMATION */
