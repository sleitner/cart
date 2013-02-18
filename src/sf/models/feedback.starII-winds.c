#include "config.h"
#ifdef STAR_FORMATION
#ifdef STAR_PARTICLE_TYPES

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

double starII_wind_momentum_boost = 1;
void starII_wind_config_init()
{
    control_parameter_add2(control_parameter_double,&starII_wind_momentum_boost,"starII:wind-momentum-boost","starII_wind_momentum_boost","factor multiplying starII_wind momentum .");
}
void starII_wind_config_verify()
{
    VERIFY(starII:wind-momentum-boost, starII_wind_momentum_boost>=0);
}

void starII_stellar_wind_kick(int level, int icell, int ipart, double ini_mass_sol, double age_yr, double Zsol, double t_next){
    double dp;
    double dt = t_next - particle_t[ipart];
    if(starII_wind_momentum_boost > 0){
	cart_assert(star_particle_type[ipart] == STAR_TYPE_STARII);
	dp = starII_wind_momentum_boost * OneStar_wind_pdot(ini_mass_sol,age_yr,Zsol) * dt;
	distribute_momentum(dp, level, icell, dt);
    }
}

#endif /* STAR_FORMATION */
#endif /* STAR_PARTICLE_TYPES */
