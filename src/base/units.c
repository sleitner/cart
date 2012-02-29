#include "config.h"

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "times.h"
#include "units.h"


struct Constants constants_internal;
const struct Constants *constants = NULL;

struct PrimaryUnits unit_factors = { 0, 0.0, 0.0, 0.0 };
const struct PrimaryUnits *primary_units = NULL;

struct Units units_internal;
const struct Units *units = NULL;


double box_size = 0.0;


#ifdef COSMOLOGY
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
#else  /* COSMOLOGY */
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


void cosmology_set_fixed();


void units_init()
{
  constants = &constants_internal;
  primary_units = &unit_factors;
  units = &units_internal;

  /*
  // Base units
  */
  constants_internal.cms = 1;
  constants_internal.erg = 1;
  constants_internal.barye = 1;
  constants_internal.dyne = 1;
  constants_internal.gpercc = 1;
  constants_internal.cc = 1;

  /* 
  //  NG: values for pc and GMsun are from http://ssd.jpl.nasa.gov/?constants
  //      values for mp, G, k, c, eV, amu are from Particle Physics Booklet 2008
  */
  constants_internal.yr = 365.25*86400; /* Julian year in seconds */
  constants_internal.Myr = 1.0e6*constants_internal.yr;
  constants_internal.Gyr = 1.0e9*constants_internal.yr;

  constants_internal.pc = 3.0856775813e18;
  constants_internal.kpc = 1.0e3*constants_internal.pc;
  constants_internal.Mpc = 1.0e6*constants_internal.pc;

  constants_internal.kms = 1.0e5;  /* cm/s */

  constants_internal.mp = 1.672621637e-24;
  constants_internal.k = 1.3806504e-16;
  constants_internal.G = 6.67428e-8;
  constants_internal.c = 2.99792458e10;

  constants_internal.eV = 1.602176487e-12;
  constants_internal.amu = 1.660538782e-24;
  constants_internal.mH  = 1.007825*constants_internal.amu;
  constants_internal.mHe = 4.002602*constants_internal.amu;

  constants_internal.Msun = 1.32712440018e26/constants_internal.G;
  constants_internal.Zsun = 0.0199; /* reference (MW) metallicity */

  constants_internal.Yp    = 0.24;                /* He mass fraction */
  constants_internal.wmu   = 4.0/(8.0-5.0*constants_internal.Yp); /* mol weight */
  constants_internal.wmu_e = 1.0/(1.0-0.5*constants_internal.Yp);
  constants_internal.XH    = 1.0 - constants_internal.Yp;
  constants_internal.XHe   = 0.25*constants_internal.Yp;
  constants_internal.gamma = 5.0/3.0;

  constants_internal.sigmaT = 6.6524e-25;

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
}


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
}
