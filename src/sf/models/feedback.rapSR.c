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

#include "starformation_feedback.h"

#include "feedback.kinetic.h"
/*
//  RaP short range from clusters
*/

float tau_UV(int icell);
static double tauUV_factor;
double rapSR_boost = 1;
double rapSR_timescale = 3e6;

double rapSR_specific_luminosity_phys;
double rapSR_timescale_phys;
double rapSR_specific_luminosity_code;
double rapSR_specific_pdot_code = 0;
double rapSR_timescale_code;
void rapSR_config_init()
{
    control_parameter_add2(control_parameter_double,&rapSR_timescale,"rapSR:timescale","rapSR_timescale","characteristic timescale over which RaP shortrange is applied.");
    control_parameter_add2(control_parameter_double,&rapSR_boost,"rapSR:boost","rapSR_boost","RaP shortrange luminosity in units of (3e42ergs/s per 1e6Msun cluster).");
}
void rapSR_config_verify()
{
    VERIFY(rapSR:boost, rapSR_boost>=0);
    VERIFY(rapSR:timescale, rapSR_timescale>0);
}
void rapSR_init()
{
    rapSR_specific_luminosity_phys = 6e42*constants->erg / (1e6*constants->Msun) ; /* rough from sb99 ergs/s per gram */
    rapSR_timescale_phys = rapSR_timescale*constants->yr;
}
void rapSR_setup(int level)
{
    tauUV_factor = units->number_density * 2.0e-21 *units->length;
    rapSR_specific_pdot_code = rapSR_specific_luminosity_phys / constants->c 
	* units->mass*units->time/(units->mass*units->velocity); /* p per time per mass */
    rapSR_specific_luminosity_code = rapSR_specific_luminosity_phys 
	/ units->energy * units->mass * units->time; /* per unit mass note this is not rt units for luminosity */
    rapSR_timescale_code = rapSR_timescale_phys / units->time ;
}

float tau_UV(int icell){
    double tau;
    tau = tauUV_factor * 
#ifdef RADIATIVE_TRANSFER
	(cell_HI_density(icell)+2.0*cell_H2_density(icell)) * rtDmw(icell) *
#else
	cell_gas_metal_density(icell) / (constants->Zsun) *
#endif
	cell_sobolev_length(icell);
    return tau;
}

double rapSR_pdot(int ipart){
#ifdef STAR_PARTICLE_TYPES
    if(!(star_particle_type[ipart] == STAR_TYPE_NORMAL 
	 || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH ) ){
	cart_error("bad particle type for rapSR");
    }
    if(rapSR_specific_pdot_code <= 0){ cart_error("rapSR not setup");}
#endif
    double tage = particle_t[ipart] - star_tbirth[ipart];
    double pdot; 
    pdot = rapSR_specific_pdot_code * particle_mass[ipart];
    if(tage > rapSR_timescale_code ){
	pdot *= pow((tage/rapSR_timescale_code),-1.2); /* fit to sb99 */
    }
    return pdot;
}

void rapSR_kick(int level, int icell, int ipart, double t_next){
    double dp, pdot;
    double dt = t_next - particle_t[ipart];
    if(rapSR_boost > 0){
	dp = rapSR_boost * rapSR_pdot(ipart) 
	    * dt * ( 1 - exp(-tau_UV(icell)) );
	distribute_momentum(dp, level, icell, dt); 
    }
}

#endif /* STAR_FORMATION */
