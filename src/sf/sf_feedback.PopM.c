#include "config.h"
#ifdef STARFORM

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
#include "starformation_feedback.h"
#include "times.h"
#include "tree.h"
#include "units.h"


extern double dUfact;
extern double feedback_temperature_ceiling;

float star_returned_advected_species[num_chem_species+1];


struct
{
  double energy_per_explosion;     /* used to be called E_51 */
  double time_duration;            /* used to be called t_fb */
  double time_delay;               
  double yield_factor;             /* fraction yield relative to the one coded in */
}
snII = { 2.0, 1.0e3, 0.0, 1.0 };


struct
{
  double energy_per_explosion;          /* used to be called E_51 */
  double time_duration;                 /* used to be called t_SNIa */
  double exploding_fraction;            /* used to be called C_SNIa */
  double mass_in_metals_per_supernova;  /* used to be called ejM_SNIa */
}
snIa = { 2.0, 2.0e8, 1.5e-2, 1.3 };


struct
{
  double loss_rate;       /* used to be called c0_ml */
  double time_interval;   /* used to be called T0_ml */
}
ml = { 0.05, 5.0e6 };


struct
{
  double ion_time_scale;
}
lum = { 3.0e6 };


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


void control_parameter_set_t0ml(const char *value, void *ptr, int ind)
{
  control_parameter_set_time(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) ml.time_interval *= 1.0e6;
}


void PopM_config_init()
{
  ControlParameterOps control_parameter_tSNIa = { control_parameter_set_tSNIa, control_parameter_list_time };
  ControlParameterOps control_parameter_t0ml = { control_parameter_set_t0ml, control_parameter_list_time };

  /*
  //  Type II supernova feedback
  */
  control_parameter_add3(control_parameter_double,&snII.energy_per_explosion,"snII:energy-per-explosion","snII.energy_per_explosion","e_51","average energy per type II supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_time,&snII.time_duration,"snII:time-duration","snII.time_duration","t_fb","time duration over which type II supernova explosions occur.");

  control_parameter_add4(control_parameter_time,&snII.time_delay,"snII:time-delay","snII.time_delay","snII:time-before-start","snII.time_before_start","time delay between the formation of the stellar particle and type II supernova explosions.");

  control_parameter_add2(control_parameter_double,&snII.yield_factor,"snII:yield-factor","snII.yield_factor","fractional yield relative to the coded in model.");

  /*
  //  Type Ia supernova feedback
  */
  control_parameter_add3(control_parameter_double,&snIa.energy_per_explosion,"snIa:energy-per-explosion","snIa.energy_per_explosion","e_51","average energy per type Ia supernova explosion, in 1e51 ergs.");

  control_parameter_add3(control_parameter_tSNIa,&snIa.time_duration,"snIa:time-duration","snIa.time_duration","t_snia","average time duration between the formation of the stellar particle and type Ia supernova explosions.");

  control_parameter_add3(control_parameter_double,&snIa.exploding_fraction,"snIa:exploding-fraction","snIa.exploding_fraction","c_snia","fraction of stars exploding as type Ia supernovae.");

  control_parameter_add3(control_parameter_double,&snIa.mass_in_metals_per_supernova,"snIa:mass-in-metals-per-supernova","snIa.mass_in_metals_per_supernova","ejm_snia","average mass (in solar masses) in metals ejected per type Ia supernova explosion.");

  /*
  //  mass loss
  */
  control_parameter_add3(control_parameter_double,&ml.loss_rate,"ml:loss-rate","ml.loss_rate","c0_ml","rate of stellar mass ejected back into the ISM by stellar mass loss. Set this to zero to disable the stellar mass loss.");

  control_parameter_add3(control_parameter_t0ml,&ml.time_interval,"ml:time-interval","ml.time_interval","t0_ml","characteristic time (in yrs) over which the stellar mass loss occurs.");

  /*
  //  ionizing luminosity
  */
  control_parameter_add2(control_parameter_time,&lum.ion_time_scale,"lum:ion-time-scale","lum.ion_time_scale","time-scale for the evolution of the ionizing luminosity from young stars.");
}


void PopM_config_verify()
{
  /*
  //  type II supernova feedback
  */
  cart_assert(!(snII.energy_per_explosion < 0.0));

  cart_assert(snII.time_duration > 0.0);

  cart_assert(!(snII.time_delay < 0.0));

  cart_assert(snII.yield_factor > 0.0);

  /*
  //  type Ia supernova feedback
  */
  cart_assert(!(snIa.energy_per_explosion < 0.0));

  cart_assert(snIa.time_duration > 0.0);

  cart_assert(snIa.exploding_fraction>0.0 && snIa.exploding_fraction<1.0);

  cart_assert(snIa.mass_in_metals_per_supernova>0.0 && snIa.mass_in_metals_per_supernova<imf->max_SNIa_mass);

  /*
  //  mass loss
  */
  cart_assert(ml.loss_rate>=0.0 && ml.loss_rate<1.0);

  cart_assert(ml.time_interval > 0.0);

  /*
  //  ionizing luminosity
  */
  cart_assert(lum.ion_time_scale > 0.0);
}


struct
{
  double energy;
  double metals;
  double teject;
  double tdelay;
}
snII_phys, snII_code;


struct
{
  double energy;
  double metals;
  double teject;
}
snIa_phys, snIa_code;


struct
{
  double dt;             /* used to be called T0_ml_code */
}
ml_code;


struct
{
  double ion_rate;
}
lum_code;


void PopM_init()
{
  int i;
  double total_mass;
  double number_SNII, number_SNIa;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( imf->fm, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_SNII = integrate( imf->f, imf->min_SNII_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(number_SNII > 0.0);

  snII_phys.tdelay = snII.time_delay;
  snII_phys.teject = snII.time_duration;
  snII_phys.energy = 1e51*constants->erg*snII.energy_per_explosion*number_SNII/(constants->Msun*total_mass);
  snII_phys.metals = snII.yield_factor*integrate( imf->fmz, imf->min_SNII_mass, imf->max_mass, 1e-6, 1e-9 )/total_mass;

  number_SNIa = snIa.exploding_fraction*integrate( imf->f, imf->min_SNIa_mass, imf->max_SNIa_mass, 1e-6, 1e-9 );
  cart_assert(number_SNIa > 0.0);

  snIa_phys.teject = snIa.time_duration;
  snIa_phys.energy = 1e51*constants->erg*snIa.energy_per_explosion*number_SNIa/(constants->Msun*total_mass);
  snIa_phys.metals = snIa.mass_in_metals_per_supernova*number_SNIa/total_mass;

  /*
  //  Decide what values of advected species stars return.
  //  Introduce a check to make sure none are forgotten.
  */
  for(i=0; i<num_chem_species; i++) star_returned_advected_species[i] = -1.1e35;

#ifdef RADIATIVE_TRANSFER
  star_returned_advected_species[(RT_HVAR_OFFSET-HVAR_ADVECTED_VARIABLES)+0] = 0.0;
  star_returned_advected_species[(RT_HVAR_OFFSET-HVAR_ADVECTED_VARIABLES)+1] = constants->XH;  /* Returned gas is ionized */
  star_returned_advected_species[(RT_HVAR_OFFSET-HVAR_ADVECTED_VARIABLES)+2] = 0.0;
  star_returned_advected_species[(RT_HVAR_OFFSET-HVAR_ADVECTED_VARIABLES)+3] = 0.0;
  star_returned_advected_species[(RT_HVAR_OFFSET-HVAR_ADVECTED_VARIABLES)+4] = constants->XHe;  /* He is doubly ionized, see Draine 2011 */
  star_returned_advected_species[(RT_HVAR_OFFSET-HVAR_ADVECTED_VARIABLES)+5] = 0.0;
#endif /* RADIATIVE_TRANSFER */

  /* Metals are returned separately, so set these to zero so as not to double count */
#ifdef ENRICHMENT
  star_returned_advected_species[HVAR_METAL_DENSITY_II-HVAR_ADVECTED_VARIABLES] = 0.0;
#ifdef ENRICHMENT_SNIa
  star_returned_advected_species[HVAR_METAL_DENSITY_Ia-HVAR_ADVECTED_VARIABLES] = 0.0;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

  /* Blaswave time is set separately */
#ifdef BLASTWAVE_FEEDBACK
  star_returned_advected_species[HVAR_BLASTWAVE_TIME-HVAR_ADVECTED_VARIABLES] = 0.0;
#endif /* BLASTWAVE_FEEDBACK*/

  /*
  //  Check that all are accounted for.
  */
  for(i=0; i<num_chem_species; i++) if(star_returned_advected_species[i] < -1.0e35)
    {
      cart_error("Advected species #%d is not accounted for in star_returned_advected_species[] array.",i);
    }

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

      cart_debug("Number of SNIa explosions per unit mass: %le per Msun",number_SNIa/total_mass);
      cart_debug("SNIa specific energy: %le erg/g = (%le km/s)^2",snIa_phys.energy,sqrt(snIa_phys.energy)/constants->kms);

#ifdef ENRICHMENT_SNIa 
      cart_debug("SNIa metal fraction : %le",snIa_phys.metals);
#endif /* ENRICHMENT_SNIa */
    }
}


void PopM_setup(int level)
{
  snII_code.tdelay = snII_phys.tdelay*constants->yr/units->time;
  snII_code.teject = snII_phys.teject*constants->yr/units->time;
  snII_code.energy = snII_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  snII_code.metals = snII_phys.metals*cell_volume_inverse[level]; 

  snIa_code.teject = snIa_phys.teject*constants->yr/units->time;
  snIa_code.energy = snIa_phys.energy*units->mass/units->energy*cell_volume_inverse[level]; 
  snIa_code.metals = snIa_phys.metals*cell_volume_inverse[level]; 

  ml_code.dt = ml.time_interval*constants->yr/units->time;

  lum_code.ion_rate = units->time/(lum.ion_time_scale*constants->yr);
}


#if defined(HYDRO) && defined(PARTICLES)

void PopM_hydrodynamic_feedback(int level, int cell, int ipart, double t_next )
{
  double dteff, phi, dU;
  double dmloss, rhor, e_old, rhofact;
  double dt = t_next - particle_t[ipart];
  double tage = particle_t[ipart] - star_tbirth[ipart];
  float thermal_pressure;
  int i;

  /* do feedback, enrichment, etc */
  if(snII_phys.energy>0.0 || snII_phys.metals>0.0)
    {
      /* snII proceeds for fpb_snII_code.teject */
      dteff = tage - snII_code.tdelay; 
      if(dteff<snII_code.teject && dteff+dt>0) 
        {
          if(dteff+dt>0 && dteff<0)
            {
              phi = min((dteff+dt)/snII_code.teject,1.0);
            }
          else
            {
              phi = min(dt,snII_code.teject-dteff)/snII_code.teject;
            }

#ifdef ENRICHMENT
          cell_gas_metal_density_II(cell) += phi*snII_code.metals*star_initial_mass[ipart];
#endif /* ENRICHMENT */

          dU = min(phi*snII_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));

          /* limit energy release and don't allow to explode in hot bubble */
          if(units->temperature*cell_gas_temperature(cell) < feedback_temperature_ceiling)
            {
              cell_gas_energy(cell) += dU;
              cell_gas_internal_energy(cell) += dU;
              cell_gas_pressure(cell) += dU*(cell_gas_gamma(cell)-1);
            }
        }
    }

  if(snIa_phys.energy>0.0 || snIa_phys.metals>0.0)
    {
      /* snIa starts at 0.1*snIa_code.dt peaks at snIa_code.dt*/
      dteff = t_next - star_tbirth[ipart];
      if(dteff > 0.1*snIa_code.teject) 
        {
          phi = f_SNIa(snIa_code.teject/dteff)*(dt/snIa_code.teject);

#ifdef ENRICHMENT_SNIa
          cell_gas_metal_density_Ia(cell) += phi*snIa_code.metals*star_initial_mass[ipart];
#endif /* ENRICHMENT_SNIa */

          dU = min(phi*snIa_code.energy*star_initial_mass[ipart],dUfact*cell_gas_density(cell));

          /* limit energy release and don't allow to explode in hot bubble */
          if(units->temperature*cell_gas_temperature(cell) < feedback_temperature_ceiling)
            {
              cell_gas_energy(cell) += dU;
              cell_gas_internal_energy(cell) += dU;
              cell_gas_pressure(cell) += dU*(cell_gas_gamma(cell)-1);
            }
        }
    }
	
  if(ml.loss_rate > 0.0)
    {
      /* limit mass loss to 10% of star's current mass */
      dmloss = min( 0.1*particle_mass[ipart],
                    star_initial_mass[ipart]*dt*ml.loss_rate / 
                    (particle_t[ipart] - (double)star_tbirth[ipart] + ml_code.dt) );
                                        
      particle_mass[ipart] -= dmloss;

      /* convert to density for cell values */
      dmloss *= cell_volume_inverse[level];

      /* account for momentum change */
      rhor = 1.0 / cell_gas_density(cell);
      e_old = cell_gas_energy(cell) -
        0.5 * ( cell_momentum(cell,0)*cell_momentum(cell,0) +
                cell_momentum(cell,1)*cell_momentum(cell,1) +
                cell_momentum(cell,2)*cell_momentum(cell,2) ) * rhor; 
      cell_gas_density(cell) += dmloss;
      rhofact = rhor * cell_gas_density(cell);
  
      cell_momentum(cell,0) += dmloss * particle_v[ipart][0];
      cell_momentum(cell,1) += dmloss * particle_v[ipart][1];
      cell_momentum(cell,2) += dmloss * particle_v[ipart][2];
                        
      cell_gas_energy(cell) = e_old + 
        0.5 * ( cell_momentum(cell,0)*cell_momentum(cell,0) +
                cell_momentum(cell,1)*cell_momentum(cell,1) +
                cell_momentum(cell,2)*cell_momentum(cell,2) ) /
        cell_gas_density(cell);

      /*
      // NG: this is to allow non-thermal pressure contribution
      */
      thermal_pressure = max((cell_gas_gamma(cell)-1.0)*cell_gas_internal_energy(cell),0.0);
      cell_gas_pressure(cell) = max(0.0,cell_gas_pressure(cell)-thermal_pressure);

      cell_gas_internal_energy(cell) *= rhofact;
      cell_gas_pressure(cell) += thermal_pressure*rhofact;

#ifdef ENRICHMENT
      cell_gas_metal_density_II(cell) += dmloss*star_metallicity_II[ipart];
#ifdef ENRICHMENT_SNIa
      cell_gas_metal_density_Ia(cell) += dmloss*star_metallicity_Ia[ipart];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
      for(i=0; i<num_chem_species; i++)
        {
          cell_advected_variable(cell,i) += dmloss*star_returned_advected_species[i];
        }
    }
}
#endif /* HYDRO && PARTICLES */


float PopM_ionizing_luminosity(int ipart)
{
  float x1, x2, dx, q, Z;

  if(!particle_is_star(ipart)) return 0.0;

  /*
  //  The convention is different from HART
  */
  x1 = lum_code.ion_rate*(particle_t[ipart]-star_tbirth[ipart]);
  if(x1 < 0.0) x1 = 0.0;
  if(x1 < 1.0e4)
    {
      /*
      //  This is a rough fit to Starburst99 evolving spectra
      */
      dx = lum_code.ion_rate*particle_dt[ipart];
      if(dx > 1.0e-5)
        {
          x2 = x1 + dx;
          x1 *= (0.8+x1*x1);
          x2 *= (0.8+x2*x2);
          q = (x2-x1)/(1+x1)/(1+x2)/particle_dt[ipart];
        }
      else
        {
          x2 = x1*(0.8+x1*x1);
          q = (0.8+3*x1*x1)/(1+x2)/(1+x2)*lum_code.ion_rate;
        }
#ifdef ENRICHMENT
      Z = (star_metallicity_II[ipart]
#ifdef ENRICHMENT_SNIa
	   + star_metallicity_Ia[ipart]
#endif /* ENRICHMENT_SNIa */
	   )/constants->Zsun;
      Z = max(1.0e-10,Z);
#else  /* ENRICHMENT */
      Z = 0.1;
#endif /* ENRICHMENT */
      return 1.04e-4/powf(Z,0.1f)/(1.0+0.27*Z)*q;     
    }
  else
    {
      return 0.0;
    }
}


struct StellarFeedback sf_feedback_PopM = 
  {
    "PopM",
    PopM_ionizing_luminosity,
    PopM_hydrodynamic_feedback,
    PopM_config_init,
    PopM_config_verify,
    PopM_init,
    PopM_setup
  };

#endif /* STARFORM */
