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


/*
//  Type Ia supernova feedback
*/

extern double dUfact;
extern double feedback_temperature_ceiling;
#ifdef ISOTROPIC_TURBULENCE_ENERGY
extern double feedback_turbulence_temperature_ceiling;
extern double fraction_SN_to_isotropic_turbulence;
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

struct
{
  double energy_per_explosion;          /* used to be called E_51 */
  double time_duration;                 /* used to be called t_SNIa */
  double exploding_fraction;            /* used to be called C_SNIa */
  double mass_in_metals_per_supernova;  /* used to be called ejM_SNIa */
  double min_mass;                      /* used to be called aM_SNIa */
  double max_mass;                 
}
  snIa = { 2.0, 2.0e8, 1.5e-2, 1.3, 3.0, 8.0 };



/* Normalized */
double f_SNIa( double x )
{
  return exp(-x*x) * sqrt(x*x*x) / 1.812804954;
}


void control_parameter_set_tSNIa(const char *value, void *ptr, int ind)
{
  control_parameter_set_time(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) snIa.time_duration *= 1.0e9;
}


void snIa_config_init()
{
  ControlParameterOps control_parameter_tSNIa = { control_parameter_set_tSNIa, control_parameter_list_time };

  control_parameter_add3(control_parameter_double,&snIa.energy_per_explosion,"snIa:energy-per-explosion","snIa.energy_per_explosion","e_51","average energy per type Ia supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_tSNIa,&snIa.time_duration,"snIa:time-duration","snIa.time_duration","t_snia","average time duration between the formation of the stellar particle and type Ia supernova explosions.");

  control_parameter_add3(control_parameter_double,&snIa.exploding_fraction,"snIa:exploding-fraction","snIa.exploding_fraction","c_snia","fraction of stars exploding as type Ia supernovae.");

  control_parameter_add3(control_parameter_double,&snIa.mass_in_metals_per_supernova,"snIa:mass-in-metals-per-supernova","snIa.mass_in_metals_per_supernova","ejm_snia","average mass (in solar masses) in metals ejected per type Ia supernova explosion.");

  control_parameter_add4(control_parameter_double,&snIa.min_mass,"snIa:min-mass","snIa.min_mass","imf:min-SNIa-mass","am_snia1","the minimum mass of stars that explode as type Ia supernovae.");

  control_parameter_add4(control_parameter_double,&snIa.max_mass,"snIa:max-mass","snIa.max_mass","imf:max-SNIa-mass","am_snia2","the maximum mass of stars that explode as type Ia supernovae.");
}


void snIa_config_verify()
{
  /*
  //  type Ia supernova feedback
  */
  VERIFY(snIa:energy-per-explosion, !(snIa.energy_per_explosion < 0.0) );

  VERIFY(snIa:time-duration, snIa.time_duration > 0.0 );

  VERIFY(snIa:exploding-fraction, snIa.exploding_fraction>0.0 && snIa.exploding_fraction<1.0 );

  VERIFY(snIa:mass-in-metals-per-supernova, snIa.mass_in_metals_per_supernova>0.0 && snIa.mass_in_metals_per_supernova<snIa.max_mass );

  VERIFY(snIa:min-mass, snIa.min_mass > 1.0 );

  VERIFY(snIa:max-mass, snIa.max_mass > snIa.min_mass );
}


struct
{
  double energy;
  double metals;
  double teject;
}
snIa_phys, snIa_code;


void snIa_init()
{
  double total_mass;
  double number_SNIa;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( imf->fm, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_SNIa = snIa.exploding_fraction*integrate( imf->f, snIa.min_mass, snIa.max_mass, 1e-6, 1e-9 );
  cart_assert(number_SNIa > 0.0);

  snIa_phys.teject = snIa.time_duration;
  snIa_phys.energy = 1e51*constants->erg*snIa.energy_per_explosion*number_SNIa/(constants->Msun*total_mass);
  snIa_phys.metals = snIa.mass_in_metals_per_supernova*number_SNIa/total_mass;

  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Number of SNIa explosions per unit mass: %le per Msun",number_SNIa/total_mass);
      cart_debug("SNIa specific energy: %le erg/g = (%le km/s)^2",snIa_phys.energy,sqrt(snIa_phys.energy)/constants->kms);

#ifdef ENRICHMENT_SNIa 
      cart_debug("SNIa metal fraction : %le",snIa_phys.metals);
#endif /* ENRICHMENT_SNIa */
    }
}


void snIa_setup(int level)
{
  snIa_code.teject = snIa_phys.teject*constants->yr/units->time;
  snIa_code.energy = snIa_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  snIa_code.metals = snIa_phys.metals*cell_volume_inverse[level]; 
}


#if defined(HYDRO) && defined(PARTICLES)

void snIa_thermal_feedback(int level, int cell, int ipart, double t_next )
{
  double dteff, phi, dU, dU_turb;

#ifdef COSMOLOGY
  double tn = tphys_from_tcode(t_next);
  double tb = tphys_from_tcode(star_tbirth[ipart]);
  double t = tphys_from_tcode(particle_t[ipart]);
#else  /* COSMOLOGY */
  double tn = t_next;
  double tb = star_tbirth[ipart];
  double t = particle_t[ipart];
#endif /* COSMOLOGY */

  double dt = tn - t;

  if(snIa_phys.energy>0.0 || snIa_phys.metals>0.0)
    {
      /* snIa starts at 0.1*snIa_code.dt peaks at snIa_code.dt*/
      dteff = tn - tb;
      if(dteff > 0.1*snIa_phys.teject) 
        {
          phi = f_SNIa(snIa_phys.teject/dteff)*(dt/snIa_phys.teject);

#ifdef ENRICHMENT_SNIa
          cell_gas_metal_density_Ia(cell) += phi*snIa_code.metals*star_initial_mass[ipart];
#endif /* ENRICHMENT_SNIa */

          dU = MIN(phi*snIa_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));
#ifdef ISOTROPIC_TURBULENCE_ENERGY
	  dU_turb = fraction_SN_to_isotropic_turbulence*dU;
	  dU_turb = MIN(  dU_turb, 
			MAX( 0.0, 
			     (feedback_turbulence_temperature_ceiling / units->temperature)
			     * cell_gas_density(cell)
			     / ((isotropic_turbulence_gamma-1)*constants->wmu) 
			     - cell_isotropic_turbulence_energy(cell)
			    )  );
	      
	  cell_isotropic_turbulence_energy(cell) += dU_turb;
	  cell_gas_energy(cell) += dU_turb;
	  cell_gas_pressure(cell) += dU_turb*(isotropic_turbulence_gamma-1);
	  dU -= dU_turb; 
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

          /* limit energy release and don't allow to explode in hot bubble */
          if(units->temperature*cell_gas_temperature(cell) < feedback_temperature_ceiling)
            {
              cell_gas_energy(cell) += dU;
              cell_gas_internal_energy(cell) += dU;
              cell_gas_pressure(cell) += dU*(cell_gas_gamma(cell)-1);
            }
        }
    }
}
#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
