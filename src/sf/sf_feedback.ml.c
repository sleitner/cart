#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "imf.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"


/*
//  Mass loss
*/
struct
{
  double loss_rate;       /* used to be called c0_ml */
  double time_interval;   /* used to be called T0_ml */
}
sfb_ml_data = { 0.05, 5.0e6 };


void control_parameter_set_t0ml(const char *value, void *ptr, int ind)
{
  control_parameter_set_time(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) sfb_ml_data.time_interval *= 1.0e6;
}


void sfb_ml_config_init()
{
  ControlParameterOps control_parameter_t0ml = { control_parameter_set_t0ml, control_parameter_list_time };

  control_parameter_add3(control_parameter_double,&sfb_ml_data.loss_rate,"ml:loss-rate","sfb_ml_data.loss_rate","c0_ml","rate of stellar mass ejected back into the ISM by stellar mass loss. Set this to zero to disable the stellar mass loss.");

  control_parameter_add3(control_parameter_t0ml,&sfb_ml_data.time_interval,"ml:time-interval","sfb_ml_data.time_interval","t0_ml","characteristic time (in yrs) over which the stellar mass loss occurs.");
}


void sfb_ml_config_verify()
{
  cart_assert(sfb_ml_data.loss_rate>=0.0 && sfb_ml_data.loss_rate<1.0);

  cart_assert(sfb_ml_data.time_interval > 0.0);
}


#if defined(HYDRO) && defined(PARTICLES)

double dt_ml_code;   /* used to be called T0_ml_code */
float star_returned_advected_species[num_chem_species+1];


void sfb_ml_init_feedback()
{
  int i;

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
}


void sfb_ml_feedback(int level, int cell, int ipart, double t_next)
{
  double dmloss, rhor, e_old, rhofact;
  double dt = t_next - particle_t[ipart];
  float thermal_pressure;
  int i;

  if(sfb_ml_data.loss_rate > 0.0)
    {
      /* limit mass loss to 10% of star's current mass */
      dmloss = min( 0.1*particle_mass[ipart],
		    star_initial_mass[ipart]*dt*sfb_ml_data.loss_rate / 
		    (particle_t[ipart] - (double)star_tbirth[ipart] + dt_ml_code) );
					
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


void sfb_ml_level_setup(int level)
{
  dt_ml_code = sfb_ml_data.time_interval*constants->yr/units->time;
}


struct StarFormationFeedback sf_feedback_ml = 
{
  "mass-loss",
  sfb_ml_feedback,
  sfb_ml_config_init,
  sfb_ml_config_verify,
  sfb_ml_init_feedback,
  sfb_ml_level_setup
};

#endif /* HYDRO && PARTICLES */
#endif /* STAR_FORMATION */
