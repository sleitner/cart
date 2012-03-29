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
//  Type II supernova feedback
*/
struct
{
  double energy_per_explosion;     /* used to be called E_51 */
  double time_duration;            /* used to be called t_fb */
  double time_delay;               
  double yield_factor;             /* fraction yield relative to the one coded in */
}
sfb_snII_data = { 2.0, 1.0e3, 0.0, 1.0 };


void sfb_snII_config_init()
{
  control_parameter_add3(control_parameter_double,&sfb_snII_data.energy_per_explosion,"snII:energy-per-explosion","sfb_snII_data.energy_per_explosion","e_51","average energy per type II supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_time,&sfb_snII_data.time_duration,"snII:time-duration","sfb_snII_data.time_duration","t_fb","time duration over which type II supernova explosions occur.");

  control_parameter_add4(control_parameter_time,&sfb_snII_data.time_delay,"snII:time-delay","sfb_snII_data.time_delay","snII:time-before-start","snII.time_before_start","time delay between the formation of the stellar particle and type II supernova explosions.");

  control_parameter_add2(control_parameter_double,&sfb_snII_data.yield_factor,"snII:yield-factor","sfb_snII_data.yield_factor","fractional yield relative to the coded in model.");
}


void sfb_snII_config_verify()
{
  cart_assert(!(sfb_snII_data.energy_per_explosion < 0.0));

  cart_assert(sfb_snII_data.time_duration > 0.0);

  cart_assert(!(sfb_snII_data.time_delay < 0.0));

  cart_assert(sfb_snII_data.yield_factor > 0.0);
}


struct
{
  double energy;
  double metals;
  double teject;
  double tdelay;
}
sfb_snII_phys, sfb_snII_code;


void sfb_snII_init_feedback()
{
  double total_mass;
  double number_SNII;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( fm_IMF, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_SNII = integrate( f_IMF, imf->min_SNII_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(number_SNII > 0.0);

  sfb_snII_phys.tdelay = sfb_snII_data.time_delay;
  sfb_snII_phys.teject = sfb_snII_data.time_duration;
  sfb_snII_phys.energy = 1e51*constants->erg*sfb_snII_data.energy_per_explosion*number_SNII/(constants->Msun*total_mass);
  sfb_snII_phys.metals = sfb_snII_data.yield_factor*integrate( fmz_IMF, imf->min_SNII_mass, imf->max_mass, 1e-6, 1e-9 )/total_mass;

  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Number of SNII explosions per unit mass: %le per Msun",number_SNII/total_mass);
      cart_debug("SNII specific energy: %le erg/g = (%le km/s)^2",sfb_snII_phys.energy,sqrt(sfb_snII_phys.energy)/constants->kms);

#ifdef ENRICHMENT 
      cart_debug("SNII metal fraction : %le",sfb_snII_phys.metals);
#endif /* ENRICHMENT */
    }
}


#if defined(HYDRO) && defined(PARTICLES)

void sfb_snII_feedback(int level, int cell, int ipart, double t_next)
{
  double dteff, phi, dU;
  double dt = t_next - particle_t[ipart];
  double tage = particle_t[ipart] - star_tbirth[ipart];

  /* do feedback, enrichment, etc */
  if(sfb_snII_phys.energy>0.0 || sfb_snII_phys.metals>0.0)
    {
      /* snII proceeds for fpb_snII_code.teject */
      dteff = tage - sfb_snII_code.tdelay; 
      if(dteff<sfb_snII_code.teject && dteff+dt>0) 
	{
	  if(dteff+dt>0 && dteff<0)
	    {
	      phi = min((dteff+dt)/sfb_snII_code.teject,1.0);
	    }
	  else
	    {
	      phi = min(dt,sfb_snII_code.teject-dteff)/sfb_snII_code.teject;
	    }

#ifdef ENRICHMENT
	  cell_gas_metal_density_II(cell) += phi*sfb_snII_code.metals*star_initial_mass[ipart];
#endif /* ENRICHMENT */

	  dU = min(phi*sfb_snII_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));

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


void sfb_snII_level_setup(int level)
{
  sfb_snII_code.tdelay = sfb_snII_phys.tdelay*constants->yr/units->time;
  sfb_snII_code.teject = sfb_snII_phys.teject*constants->yr/units->time;
  sfb_snII_code.energy = sfb_snII_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  sfb_snII_code.metals = sfb_snII_phys.metals*cell_volume_inverse[level]; 
}


struct StarFormationFeedback sf_feedback_snII = 
{
  "snII",
  sfb_snII_feedback,
  sfb_snII_config_init,
  sfb_snII_config_verify,
  sfb_snII_init_feedback,
  sfb_snII_level_setup
};

#endif /* HYDRO && PARTICLES */
#endif /* STAR_FORMATION */
