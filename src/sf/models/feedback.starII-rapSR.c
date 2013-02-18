#include "config.h"
#ifdef STAR_FORMATION
#ifdef STAR_PARTICLE_TYPES

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "imf.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "tree.h"
#include "units.h"

#include "feedback.ml.h"
#include "onestarfits.h"
#include "feedback.kinetic.h"

float tau_UV(int icell);
static double tauUV_factor;
int starII_rapSR_boost = 1;

void starII_rapSR_config_init()
{
    control_parameter_add2(control_parameter_double,&starII_rapSR_boost,"starII:rapSR-boost","starII_rapSR_boost","factor multiplying RaP shortrange.");
}
void starII_rapSR_config_verify()
{
    VERIFY(starII:rapSR-boost, starII_rapSR_boost>=0);
}
void starII_rapSR_setup(int level)
{
    tauUV_factor = units->number_density * 2.0e-21 *units->length;
}
void starII_rapSR_kick(int level, int icell, int ipart, double ini_mass_sol, double age_yr, double Zsol, double t_next){
    double dp;
    double L_UV;
    double tau;
    double dt = t_next - particle_t[ipart];
 
    cart_assert(star_particle_type[ipart] == STAR_TYPE_STARII);
    if(starII_rapSR_boost > 0){
	L_UV =  OneStar_Lbol_Lsun(ini_mass_sol,age_yr,Zsol) * OneStar_UV_fraction(ini_mass_sol, age_yr, Zsol) * Lsun_to_ergs ;
	tau = tau_UV(icell);
	
	dp  = starII_rapSR_boost * L_UV / constants->c * dt * units->time 
	    * ( 1 - exp(-tau) ) 
	    / (units->mass * units->velocity);  
	distribute_momentum(dp, level, icell, dt);

/* #ifdef DEBUG_SNL */
/* 	cart_debug("dp %e LUV %e 1-etau=%e=1",dp, dp*units->mass*units->velocity, L_UV, (1-exp(-tau)) ); */
/* #endif /\*  DEBUG_SNL *\/ */
    }
}
    
#endif /* STAR_FORMATION */
#endif /* STAR_PARTICLE_TYPES */
