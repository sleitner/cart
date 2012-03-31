#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_recipe.h"
#include "tree.h"
#include "units.h"


/*
//  Oldstyle HART recipe
*/
struct 
{
  double factor;              /* normalization constant; used to be called fmass */
  double slope;               /* used to be called alpha_SF */
  double efficiency;          /* used to be called eps_SF */
}
sfr = { 0.0, 1.5, 1.5 };


void sfr_config_init()
{
  control_parameter_add4(control_parameter_double,&sfr.slope,"sf:slope","sfr.slope","sf:recipe=0:slope","alpha_sf","the slope of the star formation law with gas density (see HART documentation for exact definition).");

  control_parameter_add4(control_parameter_double,&sfr.efficiency,"sf:efficiency","sfr.efficiency","sf:recipe=0:slope","eps_sf","the relative efficiency of the star formation law (see HART documentation for exact definition).");
}


void sfr_config_verify()
{
  cart_assert(sfr.slope > 0.0);
  cart_assert(sfr.efficiency > 0.0);
}


void sfr_setup(int level)
{
  sfr.factor = sfr.efficiency*units->time/(4.0e9*constants->yr)*pow(units->density*pow(constants->Mpc,3.0)/(1.0e16*constants->Msun),sfr.slope-1);
}


double sfr_rate(int cell)
{
  return sfr.factor*pow(cell_gas_density(cell),sfr.slope);
}


struct StarFormationRecipe sf_recipe_internal =
{
  "hart",
  sfr_rate,
  sfr_config_init,
  sfr_config_verify,
  sfr_setup
};


#endif /* HYDRO && STAR_FORMATION */
