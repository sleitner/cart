#include "config.h"

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "timestep.h"
#include "units.h"


/*
//  In principle, there is noting sacred in CGS
//  If everything is correct, the physical results 
//  should be independent of these numbers.
*/
struct CGS cgs_CGS = { 1.0, 1.0, 1.0, 1.0 };
struct CGS cgs_SI  = { 1.0e-2, 1.0e-3, 1.0, 1.0 };
extern const struct CGS *cgs = &cgs_CGS;

struct Constants constants_internal;
const struct Constants *constants = &constants_internal;

struct Units units_internal;
const struct Units *units = &units_internal;

#ifdef LEGACY_UNITS
struct LegacyUnits legacy_units_internal;
const struct LegacyUnits *legacy_units = &legacy_units_internal;
#endif /* LEGACY_UNITS */

struct PrimaryUnits unit_factors = { 0, 0.0, 0.0, 0.0 };
const struct PrimaryUnits *primary_units = &unit_factors;


double box_size = 0.0;


#ifndef COSMOLOGY

void units_set(double mass, double time, double length)
{
  cart_assert(mass > 0.0);
  cart_assert(time > 0.0);
  cart_assert(length > 0.0);

  unit_factors.mass = mass;
  unit_factors.time = time;
  unit_factors.length = length;

  unit_factors.set = 1;
}

#endif /* COSMOLOGY */


#ifdef COSMOLOGY
void control_parameter_set_OmegaM(const char *value, void *ptr, int ind)
{
  double v;
  control_parameter_set_double(value,&v,ind);
  cosmology_set(OmegaM,v);
}


void control_parameter_set_OmegaL(const char *value, void *ptr, int ind)
{
  double v;
  control_parameter_set_double(value,&v,ind);
  cosmology_set(OmegaL,v);
}


void control_parameter_set_OmegaB(const char *value, void *ptr, int ind)
{
  double v;
  control_parameter_set_double(value,&v,ind);
  cosmology_set(OmegaB,v);
}


void control_parameter_set_h(const char *value, void *ptr, int ind)
{
  double v;
  control_parameter_set_double(value,&v,ind);
  cosmology_set(h,v);
}
#endif /* COSMOLOGY */


void config_init_units()
{
#ifdef COSMOLOGY
  ControlParameterOps control_parameter_OmegaM =  { control_parameter_set_OmegaM,  control_parameter_list_double };
  ControlParameterOps control_parameter_OmegaL =  { control_parameter_set_OmegaL,  control_parameter_list_double };
  ControlParameterOps control_parameter_OmegaB =  { control_parameter_set_OmegaB,  control_parameter_list_double };
  ControlParameterOps control_parameter_h =       { control_parameter_set_h,       control_parameter_list_double };
#endif /* COSMOLOGY */

  /* 
  //  NG: values for pc and GMsun are from http://ssd.jpl.nasa.gov/?constants
  //      values for mp, G, k, c, eV, amu are from Particle Physics Booklet 2008
  */
  constants_internal.yr = 365.25*86400*cgs->s; /* Julian year in seconds */
  constants_internal.Myr = 1.0e6*constants_internal.yr;
  constants_internal.Gyr = 1.0e9*constants_internal.yr;

  constants_internal.pc = 3.0856775813e18*cgs->cm;
  constants_internal.kpc = 1.0e3*constants_internal.pc;
  constants_internal.Mpc = 1.0e6*constants_internal.pc;

  constants_internal.kms = 1.0e5*cgs->cm/cgs->s;  /* cm/s */

  constants_internal.mp = 1.672621637e-24*cgs->g;
  constants_internal.k = 1.3806504e-16*cgs->g*cgs->cm*cgs->cm/(cgs->s*cgs->s*cgs->K);
  constants_internal.G = 6.67428e-8*cgs->cm*cgs->cm*cgs->cm/(cgs->g*cgs->s*cgs->s);
  constants_internal.c = 2.99792458e10*cgs->cm/cgs->s;

  constants_internal.eV = 1.602176487e-12*cgs->g*cgs->cm*cgs->cm/(cgs->g*cgs->s);
  constants_internal.amu = 1.660538782e-24*cgs->g;
  constants_internal.mH  = 1.007825*constants_internal.amu;
  constants_internal.mHe = 4.002602*constants_internal.amu;

  constants_internal.Msun = 1.32712440018e26*cgs->cm*cgs->cm*cgs->cm/(cgs->s*cgs->s)/constants_internal.G;
  constants_internal.Zsun = 0.0199; /* reference (MW) metallicity */

  constants_internal.Yp    = 0.24;                /* He mass fraction */
  constants_internal.wmu   = 4.0/(8.0-5.0*constants_internal.Yp); /* mol weight */
  constants_internal.wmu_e = 1.0/(1.0-0.5*constants_internal.Yp);
  constants_internal.XH    = 1.0 - constants_internal.Yp;
  constants_internal.XHe   = 0.25*constants_internal.Yp;
  constants_internal.gamma = 5.0/3.0;

#ifdef COSMOLOGY

  /*
  //  Legacy parameters: cosmology MUST be set in init_run()
  */
  control_parameter_add2(control_parameter_OmegaM,(void *)&(cosmology->OmegaM),"%OmegaM","omega0","total matter density parameter.");

  control_parameter_add2(control_parameter_OmegaL,(void *)&(cosmology->OmegaL),"%OmegaL","omegal0","cosmological constant matter density parameter.");

  control_parameter_add2(control_parameter_OmegaB,(void *)&(cosmology->OmegaB),"%OmegaB","omegab0","baryonic matter density parameter.");

  control_parameter_add2(control_parameter_h,(void *)&(cosmology->h),"%h","hubble","hubble constant in 100 km/s/Mpc.");

#else  /* COSMOLOGY */

  control_parameter_add(control_parameter_double,&unit_factors.mass,"units:mass","a unit of mass in grams. This is a value that corresponds to the mass of 1 in code units.");

  control_parameter_add(control_parameter_double,&unit_factors.time,"units:time","a unit of time in seconds. This is a value that corresponds to the time interval of 1 in code units.");

  control_parameter_add(control_parameter_double,&unit_factors.length,"units:length","a unit of length in cm. This is a value that corresponds to the length of 1 in code units.");

#endif /* COSMOLOGY */

  control_parameter_add2(control_parameter_double,&box_size,"%box-size","Lbox","size of the computational box in CHIMPs.");

}


void cosmology_set_fixed();

void config_verify_units()
{
#ifdef COSMOLOGY
  if(control_parameter_is_set("%OmegaM") || control_parameter_is_set("%OmegaB") || control_parameter_is_set("%h"))
    {
      cosmology_set_fixed();
    }
#endif /* COSMOLOGY */
}


void units_set_art(double OmegaM, double h, double Lbox)
{
  double H0 = 100*constants->kms*h/constants->Mpc;

  cart_assert(OmegaM > 0.0);
  cart_assert(h > 0.0);
  cart_assert(Lbox > 0.0);

  unit_factors.length = constants->Mpc/h*Lbox/num_grid;
  unit_factors.mass = 3*pow(H0,2.0)*OmegaM/(8*M_PI*constants->G)*pow(unit_factors.length,3.0);
  unit_factors.time = 2/(H0*sqrt(OmegaM));

  unit_factors.set = 1;
}


void config_append_units_to_file(const char *filename);


void units_reset()
{
#ifdef COSMOLOGY
  if(!cosmology_is_set())
    {
      cart_error("Cosmological parameters have not been set.\nPlease set them inside the init_run() call.");
    }

  cosmology_init();
  units_set_art(cosmology->OmegaM,cosmology->h,box_size);

#else  /* COSMOLOGY */

  /* Should be already set */

  if(!unit_factors.set)
    {
      cart_error("Units has not been set in init_run().\nPlease set units with units_set(...) call.");
    }

#endif /* COSMOLOGY */


#ifdef LEGACY_UNITS
  /*
  //  Legacy units for backward compativility
  */
  /* H0 in s^-1 */
  legacy_units_internal.H0 = 100 * constants->kms * cosmology->h / constants->Mpc;

  /* r0 in h^-1 Mpc (L_box = 1.0) */
  legacy_units_internal.r0 = box_size / (double)(num_grid);

  /* t0 in years */
  legacy_units_internal.t0 = 2 / ( legacy_units_internal.H0 * sqrt(cosmology->OmegaM) ) / constants->yr;

  /* v0 in km/s */
  legacy_units_internal.v0 = 50 * legacy_units_internal.r0 * sqrt(cosmology->OmegaM);

  /* rho0 in Msun/Mpc^3 */
  legacy_units_internal.rho0 = ( 3 * legacy_units_internal.H0 * legacy_units_internal.H0 * cosmology->OmegaM / ( 8 * M_PI * constants->G ) ) * (constants->Mpc*constants->Mpc*constants->Mpc/constants->Msun);

  /* den0 in 1/cm^3 */
  legacy_units_internal.den0 = 1.123e-5 * cosmology->Omh2;

  /* P0 in g/cm/s^2 (rho0*v0^2) */
  legacy_units_internal.P0 = 4.697e-16 * cosmology->OmegaM * cosmology->OmegaM * legacy_units_internal.r0 * legacy_units_internal.r0 * cosmology->h * cosmology->h;

  /* T0 in K */
  legacy_units_internal.T0 = 3.03e5 * legacy_units_internal.r0 * legacy_units_internal.r0 * constants->wmu * cosmology->OmegaM;

  /* S0 in keV cm^2 (assumes gamma = 5/3) */
  legacy_units_internal.S0 = 52.077 * pow(constants->wmu,5./3.) * pow(cosmology->h,-4./3.) * pow(cosmology->OmegaM,1./3.) * legacy_units_internal.r0 * legacy_units_internal.r0;

  legacy_units_internal.AL_SD = 1.6625 * (1.0 - constants->Yp) * (1.0 - constants->Yp) * cosmology->h / sqrt(cosmology->OmegaM) / ( legacy_units_internal.r0 * legacy_units_internal.r0 );

  /* E0 in ergs */
  legacy_units_internal.E0 = 1.38e58 * pow(legacy_units_internal.r0,5.0) * cosmology->OmegaM * cosmology->OmegaM / cosmology->h;

  /* mass conversion */
  legacy_units_internal.M0 = legacy_units_internal.rho0 * pow(box_size/cosmology->h,3.0) / (double)num_root_cells;
#endif /* LEGACY_UNITS */

  config_append_units_to_file("config.log");
}


#ifdef LEGACY_UNITS
void check_legacy_unit(const char *name, double tol, double vnew, double vold)
{
  double d;

  if(tol < 1.0e-10) tol = 1.0e-10;

  d = fabs(vold/vnew-1.0);
  if(d > tol)
    {
      cart_error("Legacy units: large error in %s: old=%-lg,  new=%-lg, error=%-lg",name,vold,vnew,d);
    }
}
#endif /* LEGACY_UNITS */


void units_update(int level)
{
  double mb;

  if(!unit_factors.set)
    {
      cart_error("units_reset() must be called before the first call to units_update(...).");
    }

#ifdef COSMOLOGY

  units_internal.mass = unit_factors.mass;
  units_internal.time = unit_factors.time*pow(abox[level],2.0);
  units_internal.length = unit_factors.length*abox[level];

#else  /* COSMOLOGY */

  units_internal.mass = unit_factors.mass;
  units_internal.time = unit_factors.time;
  units_internal.length = unit_factors.length;

#endif /* COSMOLOGY */

  mb = constants->XH*constants->mH + constants->XHe*constants->mHe;
  
  units_internal.density = units_internal.mass/pow(units_internal.length,3.0);
  units_internal.velocity = units_internal.length/units_internal.time;
  units_internal.energy = units_internal.mass*pow(units_internal.velocity,2.0);
  units_internal.temperature = pow(units_internal.velocity,2.0)*mb/constants->k;
  units_internal.energy_density = units_internal.energy/pow(units_internal.length,3.0);
  units_internal.number_density = units_internal.density/mb;

#ifdef COSMOLOGY
  units_internal.length_in_chimps = unit_factors.length*cosmology->h/constants->Mpc;
#endif /* COSMOLOGY */

  /*
  //  Potential is slightly special; also, factor of 6 on the RHS is
  //  accounted for separately.
  */
  units_internal.potential = (2*M_PI*constants->G*units_internal.density)/(3*pow(units_internal.velocity,2.0))*pow(unit_factors.length,2.0);

  /*
  //  Other useful quantities
  */
#ifdef HYDRO
  units_internal.Emin = gas_temperature_floor/(units->temperature*constants->wmu*(constants->gamma-1));
#else
  units_internal.Emin = 0.0;
#endif /* HYDRO */

#ifdef LEGACY_UNITS
  if(level == min_level)
    {
      check_legacy_unit("mass",0.0,units->mass,legacy_units_internal.M0*constants->Msun);
      check_legacy_unit("time",0.0,units->time,legacy_units_internal.t0*pow(abox[level],2.0)*constants->yr);
      check_legacy_unit("length",0.0,units->length,legacy_units_internal.r0*abox[level]*constants->Mpc/cosmology->h);
      check_legacy_unit("energy",3.0e-4,units->energy,legacy_units_internal.E0/pow(abox[level],2.0)*cgs->g*pow(cgs->cm/cgs->s,2.0));
      check_legacy_unit("density",0.0,units->density,legacy_units_internal.rho0/pow(abox[level],3.0)*constants->Msun/pow(constants->Mpc,3.0));
      check_legacy_unit("entropy",3.0e-3,units->energy_density/pow(units->number_density,5.0/3.0),legacy_units_internal.S0*1.0e3*constants->eV/pow(cgs->cm,2.0)/pow(constants->wmu,5.0/3.0));
      check_legacy_unit("velocity",0.0,units->velocity,legacy_units_internal.v0/abox[level]*constants->kms);
      check_legacy_unit("potential",0.0,units->potential,1.0/abox[level]);
      /* T0 was rather inaccurate in legacy units */
      check_legacy_unit("temperature",2.0e-3,units->temperature,legacy_units_internal.T0/pow(abox[level],2.0)/constants->wmu*cgs->K);
      check_legacy_unit("energy_density",3.0e-4,units->energy_density,legacy_units_internal.P0/pow(abox[level],5.0)*cgs->g/cgs->cm/pow(cgs->s,2.0));
      check_legacy_unit("number_density",2.0e-3,units->number_density,legacy_units_internal.den0/pow(cgs->cm*abox[level],3.0));
      check_legacy_unit("pressure_floor_factor",0.0,constants->G*pow(units->density*units->length,2.0)/units->energy_density,abox[level]*1.5/M_PI);
#ifdef COSMOLOGY
      check_legacy_unit("length_in_chimps",0.0,units->length_in_chimps,legacy_units_internal.r0);
#endif /* COSMOLOGY */
    }
#endif /* LEGACY_UNITS */
}
