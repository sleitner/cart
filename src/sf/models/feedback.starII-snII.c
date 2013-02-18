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
#include "feedback.snII.h"

extern double dUfact;
extern double feedback_temperature_ceiling;
extern double feedback_turbulence_temperature_ceiling;
#ifdef TURBULENT_ENERGY
extern double fraction_SN_to_turbulence;
#endif /* TURBULENT_ENERGY */

struct
{
    double thermal_energy_per_explosion;
    double kinetic_energy_per_explosion;
}
starII = {1.0,0.0};
struct
{
  double thermal_energy;
  double kinetic_energy;
}
starII_SNe_phys, starII_SNe_code;

void starII_explosion_config_init()
{
    control_parameter_add2(control_parameter_double,&starII.thermal_energy_per_explosion,"starII:thermal-energy-per-explosion","starII.thermal_energy_per_explosion","thermal energy input during starII's supernova explosion, in 1e51 ergs.");
    control_parameter_add2(control_parameter_double,&starII.kinetic_energy_per_explosion,"starII:kinetic-energy-per-explosion","starII.kinetic_energy_per_explosion","kinetic energy input during starII's supernova explosion, in 1e51 ergs.");
}
void starII_explosion_config_verify()
{
  VERIFY(starII:thermal-energy-per-explosion, !(starII.thermal_energy_per_explosion < 0.0) );
  VERIFY(starII:kinetic-energy-per-explosion, !(starII.kinetic_energy_per_explosion < 0.0) );
}
void starII_explosion_feedback_init(){
    starII_SNe_phys.thermal_energy = 1e51*constants->erg*starII.thermal_energy_per_explosion; 
    starII_SNe_phys.kinetic_energy = 1e51*constants->erg*starII.kinetic_energy_per_explosion; 
}
void starII_explosion_setup(int level)
{
    starII_SNe_code.thermal_energy = starII_SNe_phys.thermal_energy / units->energy; /* not per volume */
    starII_SNe_code.kinetic_energy = starII_SNe_phys.kinetic_energy / units->energy; /* not per volume */
}
double starII_explosion_momentum(int ipart){
    /* 90% of the mass contributes to momentum. 10% stays in remnant*/
    if( particle_mass[ipart]*units->mass/constants->Msun  < snII.max_mass){
	return sqrt( 2*starII_SNe_code.kinetic_energy * 0.9*particle_mass[ipart] ); /* 2*0.5m^2v^2 */
    }else{
	return 0;
    }
}
void starII_explosion_mass(int level, int icell, int ipart){
    /* First put star back into cell: */
    putback_stellar_mass( particle_mass[ipart],level,icell,ipart );
    /* Now add type II feedback */
#ifdef ENRICHMENT
    cell_gas_metal_density_II(icell) += snII.yield_factor*OneStar_snII_Mejected_Ox(star_initial_mass[ipart]) * cell_volume_inverse[level] ; 
#endif /* ENRICHMENT */
}
void starII_explosion_kicks(int level, int icell, int ipart){
    double dp;
    if(starII_SNe_code.kinetic_energy >0.0){
	cart_assert(star_particle_type[ipart] == STAR_TYPE_STARII);
	/* this is already time-integrated momentum */
	dp = starII_explosion_momentum(ipart); 
	distribute_momentum(dp, level, icell, particle_dt[ipart]);
    }
}
void starII_explosion_thermal(int level, int icell, int ipart){
    float dU, dU_turb;
    /* Now add type II feedback */
    dU = MIN(starII_SNe_code.thermal_energy*cell_volume_inverse[level], dUfact*cell_gas_density(icell));

#ifdef TURBULENT_ENERGY
    dU_turb = fraction_SN_to_turbulence*dU;
    if(units->temperature*cell_turbulence_temperature(icell) < feedback_turbulence_temperature_ceiling)
	{
	    cell_turbulent_energy(icell) += dU_turb;
	    cell_gas_energy(icell) += dU_turb;
	    cell_gas_pressure(icell) += dU_turb*(turbulence_gamma-1);

	    dU = (1-fraction_SN_to_turbulence)*dU;
	}
#endif /* TURBULENT_ENERGY */

    /* limit energy release and don't allow to explode in hot bubble */
    if(units->temperature*cell_gas_temperature(icell) < feedback_temperature_ceiling)
	{
	    cell_gas_energy(icell) += dU;
	    cell_gas_internal_energy(icell) += dU;
	    cell_gas_pressure(icell) += dU*(cell_gas_gamma(icell)-1);
	}

}

#endif /* STAR_FORMATION */
#endif /* STAR_PARTICLE_TYPES */
