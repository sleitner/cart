#include "config.h"

#include <math.h>

#include "auxiliary.h"
#include "units.h"


#ifndef COSMOLOGY
void oldstyle_units_set(double OmegaM, double h, double Lbox)
{
  double length, mass, time;
  double H0 = 100*constants->kms*h/constants->Mpc;

  cart_assert(OmegaM > 0.0);
  cart_assert(h > 0.0);
  cart_assert(Lbox > 0.0);

  length = constants->Mpc/h*Lbox/num_grid;
  mass = 3*pow(H0,2.0)*OmegaM/(8*M_PI*constants->G)*pow(length,3.0);
  time = 2/(H0*sqrt(OmegaM));

  units_set(mass,time,length);

  units_init();
  units_update(min_level);
}
#endif
