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

extern double sf_sampling_timescale; /* just check sf_timescale is >2x longer */

double sf_timescale = 3.0e7;                 /* in yrs; used to be called tau_SF, did not exist in HART */
double sf_min_stellar_particle_mass = 0.0;   /* in Msun; used to be called dm_star_min */
double mstar_min;
double dt_SF;
void star_form_config_init()
{
    control_parameter_add4(control_parameter_double,&sf_min_stellar_particle_mass,"sf:min-stellar-particle-mass","sf_min_stellar_particle_mass","sf_min_stellar_particle_mass","dm_star_min","minimum mass for a newly created stellar particle, in solar masses. This value should be small enough to avoid artifically boosting the SFR in the low density gas.");

    control_parameter_add3(control_parameter_time,&sf_timescale,"sf:timescale","sf_timescale","tau_SF","the timescale for star formation. Star formation in a given cell is assumed to continue with the constant rate for that period of time.");
}


void star_form_config_verify()
{
    VERIFY(sf:timescale, sf_timescale >= 2*sf_sampling_timescale );

    VERIFY(sf:min-stellar-particle-mass, !(sf_min_stellar_particle_mass < 0.0) );
}


void star_form_setup(int level)
{
    mstar_min = sf_min_stellar_particle_mass * constants->Msun / units->mass; 
    
    dt_SF = sf_timescale * constants->yr / units->time;
}

void star_form_particles(int level, int icell, double dtl, double dt_eff, float sfr){
    int i;
    float mstar;
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
