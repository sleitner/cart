#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>
#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_recipe.h"
#include "rand.h"
#include "rt.h"
#include "tree.h"
#include "units.h"

extern int poissonRF12_starformation_indicator;
extern double sf_sampling_timescale; /* just check sfRF12_timescale is >2x longer */
double sfRF12_timescale = 3.0e7;                 /* in yrs; used to be called tau_SF, did not exist in HART */
double sfRF12_min_stellar_particle_mass = 0.0;   /* in Msun; used to be called dm_star_min */
double mstar_min;
double dt_SF;

void poissonRF12_config_init()
{
//    if(!(poissonRF12_starformation_indicator)) return;
    control_parameter_add4(control_parameter_time,&sfRF12_timescale,"sfRF12:timescale","sfRF12_timescale","sf_timescale","tau_SF","the timescale for star formation. Star formation in a given cell is assumed to continue with the constant rate for that period of time.");
    
    control_parameter_add4(control_parameter_double,&sfRF12_min_stellar_particle_mass,"sfRF12:min-stellar-particle-mass","sfRF12_min_stellar_particle_mass","sf_min_stellar_particle_mass","dm_star_min","minimum mass for a newly created stellar particle, in solar masses. This value should be small enough to avoid artifically boosting the SFR in the low density gas.");

}
void poissonRF12_config_verify()
{
//    if(!(poissonRF12_starformation_indicator)) return;

    VERIFY(sfRF12:timescale, sfRF12_timescale >= 2*sf_sampling_timescale );

    VERIFY(sfRF12:min-stellar-particle-mass, !(sfRF12_min_stellar_particle_mass < 0.0) );

}
void poissonRF12_init()
{
    if(!(poissonRF12_starformation_indicator)) return;
}
void poissonRF12_setup(int level)
{
    if(!(poissonRF12_starformation_indicator)) return;
    mstar_min = sfRF12_min_stellar_particle_mass * constants->Msun / units->mass; 
    
    dt_SF = sfRF12_timescale * constants->yr / units->time;
}

void poissonRF12_star_formation( int level, int icell, double dtl, double dt_eff, double sfr )
{
    if(!(poissonRF12_starformation_indicator)) return;
    int i;
    double mstar;
#ifdef OLDSTYLE_SF_ALGORITHM
    double P_SF, P_mass;
#endif
    
#ifdef OLDSTYLE_SF_ALGORITHM
    /* probability of forming a star is Poisson with <t> = dt_SF */
    P_SF = exp( -dt_eff / dt_SF );
#endif /* OLDSTYLE_SF_ALGORITHM */

    mstar = sfr*dt_SF*cell_volume[level];
    
#ifdef OLDSTYLE_SF_ALGORITHM
    
    if(mstar < mstar_min) P_mass = 1 - mstar/mstar_min; else P_mass = 0;
    
    /* randomly generate particle on timescale dt_SF */
    if(cart_rand() > P_SF+P_mass-P_SF*P_mass)
        mstar = MAX(mstar_min,mstar);
    else 
        mstar = 0.0;
    
#else  /* OLDSTYLE_SF_ALGORITHM */
    
    /* draw number of star formation events 0...\inf from poisson distribution */
    mstar = MAX( mstar_min, mstar );
    mstar *= (double)cart_rand_poisson( sfr*cell_volume[level]*dt_eff/mstar );
    
#endif /* OLDSTYLE_SF_ALGORITHM */
    
    if ( mstar > 0.0 ) {
        
        /* create the new star */
        create_star_particle( icell, mstar, dtl, STAR_TYPE_NORMAL );
        
#ifdef LOG_STAR_CREATION      
        log_star_creation( icell, particle_mass[ipart], FILE_RECORD);
#endif
        
    }
    
}

#endif /* HYDRO && STAR_FORMATION */
