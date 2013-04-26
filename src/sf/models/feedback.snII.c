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


#include "feedback.snII.h"
#include "feedback.kinetic.h"
/*
//  Type II supernova feedback
*/

extern double dUfact;
extern double feedback_temperature_ceiling;
#ifdef TURBULENT_ENERGY
extern double feedback_turbulence_temperature_ceiling;
extern double fraction_SN_to_turbulence;
#endif /* TURBULENT_ENERGY */

struct SNII_t snII = { 2.0, 1.0, 1.0e3, 0.0, 1.0, 8.0, 100.0 };
struct SNII_PROP_t snII_phys, snII_code;


void snII_config_init()
{
  control_parameter_add3(control_parameter_double,&snII.energy_per_explosion,"snII:energy-per-explosion","snII.energy_per_explosion","e_51","average energy per type II supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_time,&snII.time_duration,"snII:time-duration","snII.time_duration","t_fb","time duration over which type II supernova explosions occur.");

  control_parameter_add4(control_parameter_time,&snII.time_delay,"snII:time-delay","snII.time_delay","snII:time-before-start","snII.time_before_start","time delay between the formation of the stellar particle and type II supernova explosions.");

  control_parameter_add2(control_parameter_double,&snII.yield_factor,"snII:yield-factor","snII.yield_factor","fractional yield relative to the coded in model.");

  control_parameter_add4(control_parameter_double,&snII.min_mass,"snII:min-mass","snII.min_max","imf:min-SNII-mass","am_snii","the minimum mass of stars that explode as type II supernovae.");

  control_parameter_add2(control_parameter_double,&snII.max_mass,"snII:max-mass","snII.max_max","the maximum mass of stars that explode as type II supernovae.");

  control_parameter_add2(control_parameter_double,&snII.kinetic_energy_per_explosion,"snII:kinetic-energy-per-explosion","snII.kinetic_energy_per_explosion","average kinetic energy per type II supernova explosion, in 1e51 ergs.");
}


void snII_config_verify()
{
  /*
  //  type II supernova feedback
  */
  VERIFY(snII:energy-per-explosion, !(snII.energy_per_explosion < 0.0) );

  VERIFY(snII:time-duration, snII.time_duration > 0.0 );

  VERIFY(snII:time-delay, !(snII.time_delay < 0.0) );

  VERIFY(snII:yield-factor, snII.yield_factor > 0.0 );

  VERIFY(snII:min-mass, snII.min_mass > 1.0 );

  VERIFY(snII:max-mass, snII.max_mass > snII.min_mass );
  
  VERIFY(snII:kinetic-energy-per-explosion, !(snII.kinetic_energy_per_explosion < 0.0) );
}




void snII_init()
{
  double total_mass, ejecta_mass;
  double number_SNII;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( imf->fm, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  ejecta_mass = integrate( imf->fm, snII.min_mass, snII.max_mass, 1e-6, 1e-9 );
  number_SNII = integrate( imf->f, snII.min_mass, snII.max_mass, 1e-6, 1e-9 );
  cart_assert(number_SNII > 0.0);

  snII_phys.tdelay = snII.time_delay;
  snII_phys.teject = snII.time_duration;
  snII_phys.energy = 1e51*constants->erg*snII.energy_per_explosion*number_SNII/(constants->Msun*total_mass);
  snII_phys.metals = snII.yield_factor*integrate( imf->fmz, snII.min_mass, snII.max_mass, 1e-6, 1e-9 )/total_mass;

  snII_phys.kinetic_energy = 1e51*constants->erg*snII.kinetic_energy_per_explosion*number_SNII/(constants->Msun*total_mass); 
  snII_phys.fmass = ejecta_mass/total_mass;  
  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Number of SNII explosions per unit mass: %le per Msun",number_SNII/total_mass);
      cart_debug("SNII specific energy: %le erg/g = (%le km/s)^2",snII_phys.energy,sqrt(snII_phys.energy)/constants->kms);

#ifdef ENRICHMENT 
      cart_debug("SNII metal fraction : %le",snII_phys.metals);
#endif /* ENRICHMENT */
    }
}


void snII_setup(int level)
{
  snII_code.tdelay = snII_phys.tdelay*constants->yr/units->time;
  snII_code.teject = snII_phys.teject*constants->yr/units->time;
  snII_code.energy = snII_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  snII_code.metals = snII_phys.metals*cell_volume_inverse[level]; 

  snII_code.kinetic_energy = snII_phys.kinetic_energy * units->mass/units->energy*cell_volume_inverse[level]; 
  snII_code.fmass = snII_phys.fmass; 
}


#if defined(HYDRO) && defined(PARTICLES)

void snII_thermal_feedback(int level, int cell, int ipart, double t_next )
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
  double tage = t - tb; 

  /* do feedback, enrichment, etc */
  if(snII_phys.energy>0.0 || snII_phys.metals>0.0)
    {
      /* snII proceeds for fpb_snII_code.teject */
      dteff = tage - snII_phys.tdelay; 
      if(dteff<snII_phys.teject && dteff+dt>0) 
        {
          if(dteff+dt>0 && dteff<0)
            {
              phi = MIN((dteff+dt)/snII_phys.teject,1.0);
            }
          else
            {
              phi = MIN(dt,snII_phys.teject-dteff)/snII_phys.teject;
            }

#ifdef ENRICHMENT
          cell_gas_metal_density_II(cell) += phi*snII_code.metals*star_initial_mass[ipart];
#endif /* ENRICHMENT */

          dU = MIN(phi*snII_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));
#ifdef TURBULENT_ENERGY
	  dU_turb = fraction_SN_to_turbulence*dU;
	  if(units->temperature*cell_isotropic_turbulence_temperature(cell) < feedback_turbulence_temperature_ceiling)
	      {
		  cell_turbulent_energy(cell) += dU_turb;
		  cell_gas_energy(cell) += dU_turb;
		  cell_gas_pressure(cell) += dU_turb*(turbulence_gamma-1);
		  
		  dU = (1-fraction_SN_to_turbulence)*dU;
	      }
#endif /* TURBULENT_ENERGY */

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

void snII_kinetic_feedback(int level, int cell, int ipart, double t_next ){
    double dteff, phi;
    double dt = t_next - particle_t[ipart];
    double tage = particle_t[ipart] - star_tbirth[ipart];
    double dp, dE, dmass;
    
#ifdef STAR_PARTICLE_TYPES
    cart_assert(star_particle_type[ipart] == STAR_TYPE_NORMAL || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH );
#endif
    if(snII_code.kinetic_energy >0.0){
	dteff = tage - snII_code.tdelay;
	if(dteff<snII_code.teject && dteff+dt>0){

	    if(dteff+dt>0 && dteff<0)   {
		phi = MIN((dteff+dt)/snII_code.teject,1.0);
	    }else{
		phi = MIN(dt,snII_code.teject-dteff)/snII_code.teject;
	    }
	    
	    dE    = phi*star_initial_mass[ipart]* snII_code.kinetic_energy*cell_volume[level];
	    /* really less mass per SNe with time, but this is fine */
	    dmass = phi*star_initial_mass[ipart]* snII_code.fmass;
	    dp = sqrt( 2*dE * 0.9*dmass ); /* 2*0.5m^2v^2  ; 0.9 = remnant fraction*/
	    distribute_momentum(dp, level, cell, dt);
	}
    }
}

#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
