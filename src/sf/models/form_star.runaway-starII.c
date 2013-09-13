#include "config.h"
#ifdef STAR_FORMATION

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
#include "rand.h"

extern int starII_runaway_indicator;

double Ostar_frac_runaway = 0.5; 
double Bstar_frac_runaway = 0.1;
double Ostar_mass = 16;
double Bstar_mass = 2;
double sample_exponential(double tau);
double starII_runaway_mean_kick(double mass_msun);
int starII_runaway_velocity(double mass_code, double vadd[nDim]);


void starII_runaway_config_init()
{
//    if(!(starII_runaway_indicator)) return;
    control_parameter_add2(control_parameter_double,&Ostar_frac_runaway,"runaway:Ostar_frac","Ostar_frac_runaway","fraction of Ostars that are runaways");
    control_parameter_add2(control_parameter_double,&Bstar_frac_runaway,"runaway:Bstar_frac","Ostar_frac_runaway","fraction of Bstars that are runaways");
    control_parameter_add2(control_parameter_double,&Ostar_mass,"runaway:Ostar_mass","Ostar_mass","minimum mass of Ostars that are runaways");
    control_parameter_add2(control_parameter_double,&Bstar_mass,"runaway:Bstar_mass","Ostar_mass","minimum mass of Bstars that are runaways");
}
void starII_runaway_config_verify()
{
//    if(!(starII_runaway_indicator)) return;
    VERIFY(runaway:Ostar_frac, Ostar_frac_runaway >= 0.0 && Ostar_frac_runaway<=1.0 );
    VERIFY(runaway:Bstar_frac, Bstar_frac_runaway >= 0.0 && Bstar_frac_runaway<=1.0 );
    VERIFY(runaway:Ostar_mass, Ostar_mass >  8.0 );
    VERIFY(runaway:Bstar_mass, Bstar_mass >  0.0 && Bstar_mass < Ostar_mass );
}

double sample_exponential(double tau){
    return -tau*log(cart_rand());
}
double starII_runaway_mean_kick(double mass_msun){
    double f = constants->kms/units->velocity; 
    return 50*pow(mass_msun/33. ,0.33)*f; /* Stone 1991 */
}
int starII_runaway_velocity(double mass_code, double vadd[nDim]){
    double uni[nDim], vabs;
    int i;
    double mass_msun = mass_code * units->mass/constants->Msun;
    cart_assert(starII_runaway_indicator);

    if( mass_msun >= Bstar_mass ){
	if( (mass_msun >= Ostar_mass && cart_rand() < Ostar_frac_runaway) ||
	    (mass_msun < Ostar_mass && cart_rand() < Bstar_frac_runaway) ){

	    cart_rand_unit_vector(uni);
	    vabs = sample_exponential( starII_runaway_mean_kick(mass_msun) );
	    for(i=0;i<nDim;i++){
		vadd[i] = uni[i]*vabs;
	    }
	    
	    return 1;
	}
    }
	    
    return 0;
}
#endif /* STAR_FORMATION */
