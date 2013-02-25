#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_recipe.h"
#include "tree.h"
#include "units.h"

#include "models/form_star.all.h"
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

  control_parameter_add4(control_parameter_double,&sfr.efficiency,"sf:efficiency","sfr.efficiency","sf:recipe=0:efficiency","eps_sf","the relative efficiency of the star formation law (see HART documentation for exact definition).");
  star_form_config_init();
}


void sfr_config_verify()
{
  VERIFY(sf:slope, sfr.slope > 0.0 );
  VERIFY(sf:efficiency, sfr.efficiency > 0.0 );
  star_form_config_verify();
}


void sfr_setup(int level)
{
  sfr.factor = sfr.efficiency*units->time/(4.0e9*constants->yr)*pow(units->density*pow(constants->Mpc,3.0)/(1.0e16*constants->Msun),sfr.slope-1);
  star_form_setup( level ); 
}


double sfr_rate(int cell)
{
  return sfr.factor*pow(cell_gas_density(cell),sfr.slope);
}

void sfr_form_star_particles(int level, int icell, double dtl, double dt, float sfr){
    star_form_particles( level, icell, dtl, dt, sfr ); 
}

struct StarFormationRecipe sf_recipe_internal =
{
  "hart",
  sfr_rate,
  sfr_config_init,
  sfr_config_verify,
  sfr_setup, 
  sfr_form_star_particles
};


#endif /* HYDRO && STAR_FORMATION */
