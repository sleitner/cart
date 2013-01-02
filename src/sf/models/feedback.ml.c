#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "imf.h"
#include "particle.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"


/*
//  Stellar Mass Loss
*/

float star_returned_advected_species[num_chem_species+1];

struct
{
  double loss_rate;       /* used to be called c0_ml */
  double time_interval;   /* used to be called T0_ml */
}
  ml = { 0.05, 5.0e6 };


void control_parameter_set_t0ml(const char *value, void *ptr, int ind)
{
  control_parameter_set_time(value,ptr,ind);
  /*
  //  Backward compatibility
  */
  if(ind == 2) ml.time_interval *= 1.0e6;
}


void ml_config_init()
{
  ControlParameterOps control_parameter_t0ml = { control_parameter_set_t0ml, control_parameter_list_time };

  control_parameter_add3(control_parameter_double,&ml.loss_rate,"ml:loss-rate","ml.loss_rate","c0_ml","rate of stellar mass ejected back into the ISM by stellar mass loss. Set this to zero to disable the stellar mass loss.");

  control_parameter_add3(control_parameter_t0ml,&ml.time_interval,"ml:time-interval","ml.time_interval","t0_ml","characteristic time (in yrs) over which the stellar mass loss occurs.");
}


void ml_config_verify()
{
  /*
  //  mass loss
  */
  VERIFY(ml:loss-rate, ml.loss_rate>=0.0 && ml.loss_rate<1.0 );

  VERIFY(ml:time-interval, ml.time_interval > 0.0 );
}

void ml_init()
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
#ifdef DUST_EVOLUTION
  star_returned_advected_species[HVAR_DUST_DENSITY-HVAR_ADVECTED_VARIABLES] = 0.0;
#endif /* DUST_EVOLUTION */
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


void ml_setup(int level)
{
}


#if defined(HYDRO) && defined(PARTICLES)

void ml_feedback(int level, int cell, int ipart, double t_next )
{
  double dmloss, rhor, e_old, rhofact;
#ifdef COSMOLOGY
  double tn = tphys_from_tcode(t_next);
  double tb = tphys_from_tcode(star_tbirth[ipart]);
  double t = tphys_from_tcode(particle_t[ipart]);
#else  /* COSMOLOGY */
  double tn = t_next;
  double tb = star_tbirth[ipart];
  double t = particle_t[ipart];
#endif /* COSMOLOGY */

  float thermal_pressure;
  int i;

  if(ml.loss_rate > 0.0)
    {
      /* limit mass loss to 10% of star's current mass */
      dmloss = MIN( 0.1*particle_mass[ipart],
                    star_initial_mass[ipart]*ml.loss_rate * (tn-t) / (t - tb  + ml.time_interval) );

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
      thermal_pressure = MAX((cell_gas_gamma(cell)-1.0)*cell_gas_internal_energy(cell),0.0);
      cell_gas_pressure(cell) = MAX(0.0,cell_gas_pressure(cell)-thermal_pressure);

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


/*
//  Corrections introduced by Sam Leitner in Apr 2012
*/
void ml_snl2012_config_init()
{
  ml.loss_rate = -1;
  ml.time_interval = -1;

  ml_config_init();
}


void ml_snl2012_config_verify()
{
  /*
  //  if mass loss was not set in parameters then align it with the appropriate IMF.
  */
  /*    Leitner & Kravtsov 2011:  */
  if(ml.loss_rate == -1.0 ){
	  if(      strcmp("Salpeter",imf->name) == 0){
		  ml.loss_rate = 0.032 ;
	  }else if(strcmp("Miller-Scalo",imf->name) == 0){ //1979 
		  ml.loss_rate = 0.05 ; //left as old default -- 0.058 in paper 
	  }else if(strcmp("Chabrier",imf->name) == 0){ //Chabrier 2001~2003
		  ml.loss_rate = 0.046 ;
	  }else if(strcmp("Kroupa",imf->name) == 0){ //2001  
		  ml.loss_rate = 0.046 ;
	  }else{
		  cart_debug("IMF '%s' does not have associated IMF parameters.",imf->name);
		  cart_error("ART is terminating.");
	  }
  }

  if(ml.time_interval == -1.0 ){
	  if(      strcmp("Salpeter",imf->name) == 0){
		  ml.time_interval = 5.13e5 ;
	  }else if(strcmp("Miller-Scalo",imf->name) == 0){ //1979 
		  ml.time_interval = 5.0e6 ; //left as old default -- 6.04e6 in paper
	  }else if(strcmp("Chabrier",imf->name) == 0){ //Chabrier 2001~2003
		  ml.time_interval = 2.76e5 ;
	  }else if(strcmp("Kroupa",imf->name) == 0){ //2001  
		  ml.time_interval = 2.76e5 ;
	  }else{
		  cart_debug("IMF '%s' does not have associated IMF parameters.",imf->name);
		  cart_error("ART is terminating.");
	  }
  }

  ml_config_verify();
}

#endif /* STAR_FORMATION */
