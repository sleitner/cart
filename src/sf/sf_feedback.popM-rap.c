#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "starformation_feedback.h"
#include "rt.h"
#include "tree.h"
#include "units.h"

#include "frt/frt_c.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.rad.h"
#include "models/feedback.ml.h"


/*
//  This file adds a rudimentary RP implementation to popM-thermal
//  Only works under RADIATIVE_TRANSFER for now
*/
#ifdef RADIATIVE_TRANSFER

double rt_rp_amplt = 1.0;
double rt_rp_slope = 1.0;

struct rtRadiationPressureFactors
{
  float rf2Prad;
  float cd2tauUV;
  float cd2tauIR;
  float cd2sigma100;
}
rt_rp_factor;

void frtCall(getrpfactors)(frt_real *rf2Prad, frt_real *rho2abc);


void sfb_config_init()
{
  snII_config_init();
  snIa_config_init();
  ml_snl2012_config_init();

  control_parameter_add(control_parameter_double,&rt_rp_amplt,"@rt:rp-amplt","temporary control for testing.");
  control_parameter_add(control_parameter_double,&rt_rp_slope,"@rt:rp-slope","temporary control for testing.");
}


void sfb_config_verify()
{
  snII_config_verify();
  snIa_config_verify();
  ml_snl2012_config_verify();

  VERIFY(@rt:rp-amplt, 1 );
  VERIFY(@rt:rp-slope, 1 );
}


void sfb_init()
{
  snII_init();
  snIa_init();
  ml_init();
}


void sfb_setup(int level)
{
  frt_real rf2Prad, rho2abc;

  snII_setup(level);
  snIa_setup(level);
  rad_setup(level);
  ml_setup(level);

  /*
  //  Don't forget that n_xi is in comoving units
  */
  frtCall(getrpfactors)(&rf2Prad,&rho2abc);

  /*
  //  Factor 0.57 is the ratio of Lbol to nu*Lnu at 1000 A, per Oscar Agertz e-mail
  */
  rt_rp_factor.rf2Prad = 0.57*(constants->k)/units->energy_density*rf2Prad;
  rt_rp_factor.cd2tauUV = rho2abc;

  /*
  //  Units for tauIR factor, assuming kIR = 5 cm^2/g, as in Hopkins et al 1101.4940
  */
  rt_rp_factor.cd2tauIR = 5*units->density*units->length;

  rt_rp_factor.cd2sigma100 = units->density*units->length/(100*constants->Msun/constants->pc/constants->pc);
}


#if defined(HYDRO) && defined(PARTICLES)
void sfb_hydro_feedback(int level, int cell, int ipart, double t_next )
{
  snII_thermal_feedback(level,cell,ipart,t_next);
  snIa_thermal_feedback(level,cell,ipart,t_next);
  ml_feedback(level,cell,ipart,t_next);
}
#endif /* HYDRO && PARTICLES */


float sfb_radiation_pressure(int cell)
{
  float len = cell_sobolev_length(cell);
  float cd_gas = len*constants->XH*cell_gas_density(cell);
  float cd_dust = len*constants->XH*rtDmw(cell)*(cell_HI_density(cell)+2*cell_H2_density(cell));
  float tauUV = rt_rp_factor.cd2tauUV*cd_dust;
  float pressure ;
#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD==RT_METHOD_OTVET)
  float rfLoc;

#ifdef RT_UV
  float w;
  if(cell_var(cell,rt_field_offset+rt_num_freqs-1)>0.0 && cell_var(cell,rt_field_offset+rt_num_freqs-1)>1.0e-35*cell_var(cell,RT_VAR_OT_FIELD))
    {
      w = log(cell_var(cell,RT_VAR_OT_FIELD)/cell_var(cell,rt_field_offset+rt_num_freqs-1));
      if(tauUV > w) tauUV = w;
    }
  rfLoc = cell_var(cell,rt_field_offset+rt_num_freqs-1)*exp(tauUV);
#else /* RT_UV */
  rfLoc = cell_var(cell,RT_VAR_OT_FIELD);
#endif /* RT_UV */
 
  pressure = (1-exp(-tauUV)+rt_rp_amplt*pow(rt_rp_factor.cd2sigma100*cd_gas,rt_rp_slope))*rt_rp_factor.rf2Prad*rfLoc;
  PLUGIN_POINT(RadiationFeedbackEnd)(level,cell,pressure);

  return pressure;

#else /* RT_TRANSFER && RT_TRANSFER_METHOD==RT_METHOD_OTVET */
  cart_error("Radiation pressure without RT_TRANSFER and RT_TRANSFER_METHOD=RT_METHOD_OTVET is not implemented yet.");
  return 0.0;
#endif /* RT_TRANSFER && RT_TRANSFER_METHOD==RT_METHOD_OTVET */
}



struct StellarFeedbackParticle sf_feedback_particle_internal =
  {
    "popM-rap",
    sfb_hydro_feedback,
    rad_luminosity_popM,
    sfb_radiation_pressure,
    sfb_config_init,
    sfb_config_verify,
    sfb_init,
    sfb_setup
    sfb_setup,
    NULL
  };

void sfb_hydro_feedback_cell(int level, int cell, double t_next, double dt ){}

struct StellarFeedbackCell sf_feedback_cell_internal =
{
    sfb_hydro_feedback_cell,
    gather_kicks
};



#else  /* RADIATIVE_TRANSFER */

#error "This implementation requires RADIATIVE_TRANSFER to be set."

#endif /* RADIATIVE_TRANSFER */

#endif /* STAR_FORMATION */
