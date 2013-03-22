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

#include "gsl/gsl_interp.h"

/*
//  feedback from AGB stars (metals & dust)
//  R. Feldmann 2012
//  stars are assumed to go through an instantaneous AGB phase at the end of their life time
*/

struct
{
  double yield_factor;             /* fractional yield relative to the one coded in */
  double dust_yield_factor;        /* fractional dust yield relative to the one coded in */
  double min_mass;                 
  double max_mass;                 
}
AGB = { 1.0, 1.0, 1.0, 8.0 };


void AGB_config_init()
{
  control_parameter_add2(control_parameter_double,&AGB.yield_factor,"AGB:yield-factor","AGB.yield_factor","fractional yield relative to the coded in model.");

  control_parameter_add2(control_parameter_double,&AGB.dust_yield_factor,"AGB:dust-yield-factor","AGB.dust_yield_factor","fractional dust yield relative to the coded in model.");

  control_parameter_add4(control_parameter_double,&AGB.min_mass,"AGB:min-mass","AGB.min_mass","imf:min-AGB-mass","am_AGB","the minimum mass of stars that enter AGB phase.");

  control_parameter_add2(control_parameter_double,&AGB.max_mass,"AGB:max-mass","AGB.max_mass","the maximum mass of stars that enter AGB phase.");
}


void AGB_config_verify()
{
  
  VERIFY(AGB:yield-factor, AGB.yield_factor > 0.0 );

  VERIFY(AGB:dust-yield-factor, AGB.dust_yield_factor > 0.0 );

  VERIFY(AGB:min-mass, AGB.min_mass >= 1.0 );

  VERIFY(AGB:max-mass, AGB.max_mass > AGB.min_mass );
}

/* extra - auxiliary data */
struct
{ 
  double minAge;
  double maxAge;
  double conv_facM;
  double conv_facD;
}
AGB_tab;

void AGB_init()
{
  double total_mass;
  double number_AGB;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( imf->fm, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_AGB = integrate( imf->f, AGB.min_mass, AGB.max_mass, 1e-6, 1e-9 );
  cart_assert(number_AGB > 0.0);

  AGB_tab.maxAge = pow(10.,tlf(log10(AGB.min_mass),log10(0.02)));
  AGB_tab.minAge = pow(10.,tlf(log10(AGB.max_mass),log10(0.02)));
  cart_assert(AGB_tab.minAge>6.6 && AGB_tab.maxAge<1e10);
  
  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE) {
    cart_debug("AGB: minMass = %g maxMass = %g [Msun]; minAge = %g maxAge = %g [yr]",AGB.min_mass, AGB.max_mass, AGB_tab.minAge, AGB_tab.maxAge);
  }
}

void AGB_setup(int level)
{
  AGB_tab.conv_facM = AGB.yield_factor * cell_volume_inverse[level];
  AGB_tab.conv_facD = AGB.dust_yield_factor * cell_volume_inverse[level];
  if(local_proc_id == MASTER_NODE) {
    /* cart_debug("AGB: level=%d conv_facM=%g conv_facD=%g",level,AGB_tab.conv_facM,AGB_tab.conv_facD); */
  }
}

#if defined(HYDRO) && defined(PARTICLES)

void AGB_feedback(int level, int cell, int ipart, double t_next )
{

#ifdef ENRICHMENT
  
  double tn, tb, t, dt;
  double lage;
  double rmetal, rdust;
 
  tn = tphys_from_tcode(t_next);
  t = tphys_from_tcode(particle_t[ipart]);
  tb =  tphys_from_tcode(star_tbirth[ipart]);
  
  if (t-tb > AGB_tab.maxAge || tn-tb < AGB_tab.minAge) /* stellar population is too old / too young for AGB */
    return; 
  dt = tn - t;
  lage = log10(t-tb);
  
  /* the following fitting functions are described in Feldmann+ in prep */
    
  /* metal injection rate [Msun yr^-1] per solar mass SSP */
  if (lage < 7.56)
    rmetal = pow(10, -7.531 -  3*pow( (lage-6.6), 0.711 ) );
  else
    rmetal = pow(10, -10.446 - 1.33*pow( (lage-7.56), 0.94 ) );
  
  /* dust injection rate [Msun yr^-1] per solar mass SSP */
  rdust = 0.1*rmetal;   

  /*
  if (level==9 && lage<8 && lage>7.5)
    printf("AGB: ipart=%d cell=%d dt=%g log_age=%g M*=%g=%g [Msun] BEFORE: gas_den=%g metal_den=%g dust_den=%g rmetal=%g rdust=%g\n",
                ipart,cell,dt,lage, star_initial_mass[ipart], star_initial_mass[ipart]*units->mass/constants->Msun,
                cell_gas_density(cell), cell_gas_metal_density(cell), cell_dust_density(cell), rmetal, rdust);
  */

  cell_gas_metal_density_II(cell) +=  rmetal * dt * star_initial_mass[ipart] * AGB_tab.conv_facM;

  #ifdef DUST_EVOLUTION
  #ifndef DUST_EVOLUTION_noAGB

  cell_dust_density(cell) +=  rdust * dt * star_initial_mass[ipart] * AGB_tab.conv_facD;

  #endif
  #endif

  /*
  if (level==9 && lage<8 && lage>7.5)
    printf("AGB: ipart=%d cell=%d dt=%g log_age=%g M*=%g AFTER: gas_den=%g metal_den=%g dust_den=%g rmetal=%g rdust=%g\n",
                ipart,cell,dt,lage, star_initial_mass[ipart]*units->mass/constants->Msun,
                cell_gas_density(cell), cell_gas_metal_density(cell), cell_dust_density(cell), rmetal, rdust);
  */

#endif /* ENRICHMENT */ 

}
#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
