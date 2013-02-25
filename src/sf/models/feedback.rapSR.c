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

double rapSR_luminosity_phys;
double rapSR_timescale_phys;
double rapSR_luminosity_code;
double rapSR_pdot_code;
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
    rapSR_luminosity_phys = rapSR_boost * 6e42*constants->erg / (1e6*constants->Msun) ; /* rough from sb99 ergs/s per gram */
    rapSR_timescale_phys = rapSR_timescale*constants->yr;
}
void rapSR_setup(int level)
{
    tauUV_factor = units->number_density * 2.0e-21 *units->length;
    rapSR_pdot_code = rapSR_luminosity_phys / constants->c 
	* units->mass*units->time/(units->mass*units->velocity); /* p per time per mass! */
    rapSR_luminosity_code = rapSR_luminosity_phys / units->energy * units->mass * units->time; /* per unit mass note this is not rt units for luminosity */
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

void rapSR_kick(int level, int icell, int ipart, double t_next){
    double dp ;
    double tau;
    double dt = t_next - particle_t[ipart];
    double tage = particle_t[ipart] - star_tbirth[ipart];
#ifdef STAR_PARTICLE_TYPES
    cart_assert(star_particle_type[ipart] == STAR_TYPE_NORMAL || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH );
#endif
    if(rapSR_boost > 0){
	tau = tau_UV(icell);
	dp  =  rapSR_pdot_code * particle_mass[ipart] * dt * ( 1 - exp(-tau) );
	if(tage > rapSR_timescale_code ){
	    dp *= pow((tage/rapSR_timescale_code),-1.2); /* fit to sb99 */
	}
#ifdef DEBUG_SNL
	double L_UV = rapSR_luminosity_phys * particle_mass[ipart]*units->mass; //ergs/s
	double L_UVRT = sf_feedback->rt_source(ipart) * particle_mass[ipart] 
	    *units->mass*constants->c*constants->c/units->time ; // ergs/s 
	double dp2   =  L_UV/constants->c * dt*units->time * ( 1 - exp(-tau) ) 
	    /(units->mass*units->velocity) ;
	double dpRT  =  L_UVRT / constants->c * dt*units->time * ( 1 - exp(-tau) )
	    /(units->mass*units->velocity) ;
	double fact = 1.0/cell_gas_density(icell)*units->velocity/constants->kms*cell_volume_inverse[level];
	cart_debug("L_UV %e, RT L_UVRT %e ||dp %e dp2 %e  dpRT %e || mass=%e dp_code %e",
		   L_UV, L_UVRT, 
		   dp*fact,  
		   dp2*fact,
		   dpRT*fact, 
		   cell_gas_density(icell)*cell_volume[level]*units->mass/constants->Msun,
		   dp
	    );
 	cart_debug("snl: starII_popM_UV_momentum");   
#endif 
	distribute_momentum(dp, level, icell, dt); 
    }
}

#endif /* STAR_FORMATION */
