#include "config.h"

#include <stdio.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "times.h"
#include "units.h"


int drive()
{
  cosmology_set_OmegaM(0.3);
  cosmology_set_h(0.7);
  cosmology_set_OmegaB(0.04);
  box_size = 1;

  units_init();

  abox[min_level] = 1;
  units_update(min_level);

  cart_debug("Unit of density = %lg",units->density);

  return 0;
}
