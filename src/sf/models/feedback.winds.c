#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "tree.h"
#include "units.h" 

#include "onestarfits.h"
#include "feedback.kinetic.h"

double wind_momentum_boost = 1;
double wind_duration = 5e6; /* Leitherer+92 */
double wind_duration_phys;
double wind_duration_code;
void wind_config_init()
{
    control_parameter_add2(control_parameter_double,&wind_momentum_boost,"wind:momentum-boost","wind_momentum_boost","factor multiplying wind momentum .");
    control_parameter_add2(control_parameter_time,&wind_duration,"wind:time-duration","wind_duration","time duration over which stellar winds occur.");
}
void wind_config_verify()
{
    VERIFY(wind:momentum-boost, wind_momentum_boost >= 0.0);
    VERIFY(wind:time-duration, wind_duration > 0.0);
}
void wind_init()
{
    wind_duration_phys = wind_duration;
}

static double sb99_wind_pdot; 
void wind_setup(int level)
{
    /* pdot in dyne per mass Leitherer+92 figure 13*/
    sb99_wind_pdot = wind_momentum_boost * 2e32 *units->time/(units->velocity*units->mass) / (1.0e6*constants->Msun/units->mass); 
    wind_duration_code = wind_duration_phys*constants->yr/units->time;
}

void stellar_wind_kick(int level,int icell,int ipart,double t_next){
    double dteff, dt, tage;
    double dp;
#ifdef STAR_PARTICLE_TYPES
    cart_assert(star_particle_type[ipart] == STAR_TYPE_NORMAL || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH );
#endif
    if(wind_momentum_boost > 0){
	dt = t_next - particle_t[ipart];
	tage = particle_t[ipart] - star_tbirth[ipart];
	
	/* winds fill in the period before supernovae */
	if(tage < wind_duration_code){
	    dteff = MIN(dt,wind_duration_code-tage);
	    dp = sb99_wind_pdot 
#ifdef ENRICHMENT
		* star_metallicity_II[ipart]/constants->Zsun 
#endif /* ENRICHMENT */
		* dteff * star_initial_mass[ipart]  ;
	
	    distribute_momentum(dp, level,  icell, dt); 
	}
    }
}
/* const double sb_wind_lum ; */
/* const double sb_wind_fmdot ; */
/*     sb_wind_lum = 2e34*constants->erg/units->energy/constants->Msun*units->mass; /\* ergs/s/Msun figure 3 *\/ */
/*     sb_wind_fmdot = 1e-8 *units->time/constants->yr; /\* mass lost per unit mass per year fig 7b *\/ */
/* 	mwind = sb_wind_fmdot * star_particle_mass[ipart] * dteff; */
/* 	dE = phi * sb_wind_lum * star_particle_mass[ipart]  * dteff; */
/* 	dp = sqrt( 2*dE*mwind ); */

#endif /* STAR_FORMATION */
