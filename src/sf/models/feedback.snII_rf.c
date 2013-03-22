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
//  feedback from Type II supernova (metals & dust)
//  R. Feldmann 2012
*/


extern double dUfact;
extern double feedback_temperature_ceiling;

struct
{
  double energy_per_explosion;     /* used to be called E_51 */
  double yield_factor;             /* fraction yield relative to the one coded in */
  double dust_yield_factor;        /* fractional dust yield relative to the one coded in */
  double dust_destruction_factor;  /* destruction efficiency relative to the one coded in */
  double min_mass;                 /* used to be called aM_snII_rf */
  double max_mass;                 
}
  snII_rf = { 1.0, 1.0, 1.0, 1.0, 8.0, 40.0 };


void snII_rf_config_init()
{
  control_parameter_add3(control_parameter_double,&snII_rf.energy_per_explosion,"snII_rf:energy-per-explosion","snII_rf.energy_per_explosion","e_51","average energy per type II supernova explosion, in 1e51 ergs.");

  control_parameter_add2(control_parameter_double,&snII_rf.yield_factor,"snII_rf:yield-factor","snII_rf.yield_factor","fractional yield relative to the one coded in model.");
  
  control_parameter_add2(control_parameter_double,&snII_rf.dust_yield_factor,"snII_rf:dust-yield-factor","snII_rf.dust_yield_factor","fractional dust yield relative to the one coded in model.");
  
  control_parameter_add2(control_parameter_double,&snII_rf.dust_destruction_factor,"snII_rf:dust-destruction-factor","snII_rf.dust_destruction_factor","dust destruction efficiency relative to the one coded in model.");

  control_parameter_add4(control_parameter_double,&snII_rf.min_mass,"snII_rf:min-mass","snII_rf.min_max","imf:min-snII_rf-mass","am_snII_rf","the minimum mass of stars that explode as type II supernovae.");

  control_parameter_add2(control_parameter_double,&snII_rf.max_mass,"snII_rf:max-mass","snII_rf.max_max","the maximum mass of stars that explode as type II supernovae.");
}


void snII_rf_config_verify()
{
  /*
  //  type II supernova feedback
  */
  VERIFY(snII_rf:energy-per-explosion, !(snII_rf.energy_per_explosion < 0.0) );

  VERIFY(snII_rf:yield-factor, snII_rf.yield_factor > 0.0 );
  
  VERIFY(snII_rf:dust-yield-factor, snII_rf.dust_yield_factor > 0.0 );
  
  VERIFY(snII_rf:dust-destruction-factor, snII_rf.dust_destruction_factor > 0.0 );

  VERIFY(snII_rf:min-mass, snII_rf.min_mass > 1.0 );

  VERIFY(snII_rf:max-mass, snII_rf.max_mass > snII_rf.min_mass );
}

struct
{ 
  double minAge;
  double maxAge;
  double conv_facE;
  double conv_facM;
  double conv_facD;
  double conv_facD2;
}
snII_rf_tab;

void snII_rf_init()
{
  double total_mass;
  double number_snII_rf;

  /*
  //  All masses are in Msun
  */  
  total_mass = integrate( imf->fm, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
  cart_assert(total_mass > 0.0);

  number_snII_rf = integrate( imf->f, snII_rf.min_mass, snII_rf.max_mass, 1e-6, 1e-9 );
  cart_assert(number_snII_rf > 0.0);

  snII_rf_tab.maxAge = pow(10.,tlf(log10(snII_rf.min_mass),log10(0.02)));
  snII_rf_tab.minAge = pow(10.,tlf(log10(snII_rf.max_mass),log10(0.02)));
  cart_assert(snII_rf_tab.minAge>6.6 && snII_rf_tab.maxAge<1e10);
 
  /*
  // The rest is for diagnostic only
  */
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("snII_rf: minMass = %g maxMass = %g [Msun]; minAge = %g maxAge = %g [yr]",snII_rf.min_mass, snII_rf.max_mass, snII_rf_tab.minAge, snII_rf_tab.maxAge);
      cart_debug("Number of snII explosions per unit mass: %le per Msun",number_snII_rf/total_mass);

#ifdef ENRICHMENT 
#endif /* ENRICHMENT */
    }
}


void snII_rf_setup(int level)
{
  snII_rf_tab.conv_facE = snII_rf.energy_per_explosion    * 1e51 *constants->erg / units->energy * units->mass / constants->Msun * cell_volume_inverse[level];
  snII_rf_tab.conv_facM = snII_rf.yield_factor * cell_volume_inverse[level];
  snII_rf_tab.conv_facD = snII_rf.dust_yield_factor * cell_volume_inverse[level];
  snII_rf_tab.conv_facD2 = snII_rf.dust_destruction_factor * 1000 * cell_volume_inverse[level]; /* default is epsilon*mass_SN = 1000 Msun */
  if(local_proc_id == MASTER_NODE)
    {
      /* cart_debug("snII_rf: level=%d conv_facE=%g conv_facM=%g conv_facD=%g conv_facD2=%g",level,snII_rf_tab.conv_facE,snII_rf_tab.conv_facM,snII_rf_tab.conv_facD,snII_rf_tab.conv_facD2); */
    }

}

#if defined(HYDRO) && defined(PARTICLES)

void snII_rf_feedback(int level, int cell, int ipart, double t_next )
{
  double dt;
  double lage;
  double rSN, rmetal, rdust;
   
  /* times in yr */
#ifdef COSMOLOGY
  double tn = tphys_from_tcode(t_next);
  double tb = tphys_from_tcode(star_tbirth[ipart]);
  double t = tphys_from_tcode(particle_t[ipart]);
#else  /* COSMOLOGY */
  double tn = t_next;
  double tb = star_tbirth[ipart];
  double t = particle_t[ipart];
#endif /* COSMOLOGY */
  
  if (t-tb > snII_rf_tab.maxAge || tn-tb < snII_rf_tab.minAge) /* stellar population is too old / too young for SN II */
    return;  
  
  dt = tn - t;
  lage = log10(t-tb);
  
  /* the following fitting functions are described in Feldmann+ in prep */
    
  /* SN rate [yr^-1] per solar mass SSP */
  rSN =  pow(10, -9.2 - 0.4*pow( (lage-6.5), 1.25 ) );
  
  /* metal injection rate [Msun yr^-1] per solar mass SSP */
  if (lage < 7.56)
    rmetal = pow(10, -7.531 -  3*pow( (lage-6.6), 0.711 ) );
  else
    rmetal = pow(10, -10.446 - 1.33*pow( (lage-7.56), 0.94 ) );
  
  /* dust injection rate [Msun yr^-1] per solar mass SSP */
  rdust = 0.1*rmetal;   

  /* energy injection */
  if (snII_rf.energy_per_explosion > 0.0) {
    
    double dU = MIN(rSN * dt * star_initial_mass[ipart] * snII_rf_tab.conv_facE, dUfact*cell_gas_density(cell) );

    if(units->temperature*cell_gas_temperature(cell) < feedback_temperature_ceiling) {
      cell_gas_energy(cell) += dU;
      cell_gas_internal_energy(cell) += dU;
      cell_gas_pressure(cell) += dU*(cell_gas_gamma(cell)-1);
    }
  }

  /* metal injection */
#ifdef ENRICHMENT
  
  /*
  if (level==9 && lage<7.5)
    printf("SNII: ipart=%d cell=%d dt=%g log_age=%g M*=%g=%g [Msun] BEFORE: gas_den=%g metal_den=%g dust_den=%g rSN=%g rmetal=%g rdust=%g\n",
                ipart,cell,dt,lage, star_initial_mass[ipart], star_initial_mass[ipart]*units->mass/constants->Msun,
                cell_gas_density(cell), cell_gas_metal_density(cell), cell_dust_density(cell), rSN, rmetal, rdust);
  */
 
  cell_gas_metal_density_II(cell) +=  rmetal * dt * star_initial_mass[ipart] * snII_rf_tab.conv_facM;

  #ifdef DUST_EVOLUTION
  #ifndef DUST_EVOLUTION_noSNII

  cell_dust_density(cell) +=  ( rdust * snII_rf_tab.conv_facD - rSN * snII_rf_tab.conv_facD2 * cell_dust_density(cell)/cell_gas_density(cell) ) * dt * star_initial_mass[ipart];
  
  if (cell_dust_density(cell)<0)
    cell_dust_density(cell)=0.;

  #endif
  #endif
  
  /* 
  if (level==9 && lage<7.5)
    printf("SNII: ipart=%d cell=%d dt=%g log_age=%g M*=%g AFTER: gas_den=%g metal_den=%g dust_den=%g rSN=%g rmetal=%g rdust=%g\n",
                ipart,cell,dt,lage, star_initial_mass[ipart]*units->mass/constants->Msun,
                cell_gas_density(cell), cell_gas_metal_density(cell), cell_dust_density(cell), rSN, rmetal, rdust);
  */
		
#endif /* ENRICHMENT */ 
    
}
#endif /* HYDRO && PARTICLES */

#endif /* STAR_FORMATION */
