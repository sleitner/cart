#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "imf.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"


extern double dUfact;
extern double feedback_temperature_ceiling;

/*
//  Type Ia supernova feedback
*/
struct
{
  double energy_per_explosion;          /* used to be called E_51 */
  double time_duration;                    /* used to be called t_SNIa */
  double exploding_fraction;            /* used to be called C_SNIa */
  double mass_in_metals_per_supernova;  /* used to be called ejM_SNIa */
}
sfb_snIa_data = { 2.0, 2.0e8, 1.5e-2, 1.3 };


void control_parameter_set_tSNIa(const char *value, void *ptr, int ind)
{
  control_parameter_set_time(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) sfb_snIa_data.time_duration *= 1.0e9;
}


void sfb_snIa_config_init()
{
  ControlParameterOps control_parameter_tSNIa = { control_parameter_set_tSNIa, control_parameter_list_time };

  control_parameter_add3(control_parameter_double,&sfb_snIa_data.energy_per_explosion,"snIa:energy-per-explosion","sfb_snIa_data.energy_per_explosion","e_51","average energy per type Ia supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_tSNIa,&sfb_snIa_data.time_duration,"snIa:time-duration","sfb_snIa_data.time_duration","t_snia","average time duration between the formation of the stellar particle and type Ia supernova explosions.");

  control_parameter_add3(control_parameter_double,&sfb_snIa_data.exploding_fraction,"snIa:exploding-fraction","sfb_snIa_data.exploding_fraction","c_snia","fraction of stars exploding as type Ia supernovae.");

  control_parameter_add3(control_parameter_double,&sfb_snIa_data.mass_in_metals_per_supernova,"snIa:mass-in-metals-per-supernova","sfb_snIa_data.mass_in_metals_per_supernova","ejm_snia","average mass (in solar masses) in metals ejected per type Ia supernova explosion.");
}


void sfb_snIa_config_verify()
{
  cart_assert(!(sfb_snIa_data.energy_per_explosion < 0.0));

  cart_assert(sfb_snIa_data.time_duration > 0.0);

  cart_assert(sfb_snIa_data.exploding_fraction>0.0 && sfb_snIa_data.exploding_fraction<1.0);

  cart_assert(sfb_snIa_data.mass_in_metals_per_supernova>0.0 && sfb_snIa_data.mass_in_metals_per_supernova<imf->max_SNIa_mass);
}


struct
{
  double energy;
  double metals;
  double teject;
}
sfb_snIa_phys, sfb_snIa_code;


/* Normalized */
double f_SNIa( double x )
{
  return exp(-x*x) * sqrt(x*x*x) / 1.812804954;
}


void sfb_snIa_init_feedback()
{
  double total_mass;
  double number_SNIa;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( fm_IMF, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_SNIa = sfb_snIa_data.exploding_fraction*integrate( f_IMF, imf->min_SNIa_mass, imf->max_SNIa_mass, 1e-6, 1e-9 );
  cart_assert(number_SNIa > 0.0);

  sfb_snIa_phys.teject = sfb_snIa_data.time_duration;
  sfb_snIa_phys.energy = 1e51*constants->erg*sfb_snIa_data.energy_per_explosion*number_SNIa/(constants->Msun*total_mass);
  sfb_snIa_phys.metals = sfb_snIa_data.mass_in_metals_per_supernova*number_SNIa/total_mass;

  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Number of SNIa explosions per unit mass: %le per Msun",number_SNIa/total_mass);
      cart_debug("SNIa specific energy: %le erg/g = (%le km/s)^2",sfb_snIa_phys.energy,sqrt(sfb_snIa_phys.energy)/constants->kms);

#ifdef ENRICHMENT_SNIa 
      cart_debug("SNIa metal fraction : %le",sfb_snIa_phys.metals);
#endif /* ENRICHMENT_SNIa */
    }
}


#if defined(HYDRO) && defined(PARTICLES)

void sfb_snIa_feedback(int level, int cell, int ipart, double t_next)
{
  double dteff, phi, dU;
  double dt = t_next - particle_t[ipart];

  if(sfb_snIa_phys.energy>0.0 || sfb_snIa_phys.metals>0.0)
    {
      /* snIa starts at 0.1*sfb_snIa_code.dt peaks at sfb_snIa_code.dt*/
      dteff = t_next - star_tbirth[ipart];
      if(dteff > 0.1*sfb_snIa_code.teject) 
	{
	  phi = f_SNIa(sfb_snIa_code.teject/dteff)*(dt/sfb_snIa_code.teject);

#ifdef ENRICHMENT_SNIa
	  cell_gas_metal_density_Ia(cell) += phi*sfb_snIa_code.metals*star_initial_mass[ipart];
#endif /* ENRICHMENT_SNIa */

	  dU = min(phi*sfb_snIa_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));

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


void sfb_snIa_level_setup(int level)
{
  sfb_snIa_code.teject = sfb_snIa_phys.teject*constants->yr/units->time;
  sfb_snIa_code.energy = sfb_snIa_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  sfb_snIa_code.metals = sfb_snIa_phys.metals*cell_volume_inverse[level]; 
}


struct StarFormationFeedback sf_feedback_snIa = 
{
  "snIa",
  sfb_snIa_feedback,
  sfb_snIa_config_init,
  sfb_snIa_config_verify,
  sfb_snIa_init_feedback,
  sfb_snIa_level_setup
};

#endif /* HYDRO && PARTICLES */
#endif /* STAR_FORMATION */
