#include "config.h"
#if defined(HYDRO) && defined(STARFORM)

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_recipe.h"
#include "tree.h"
#include "units.h"


/*
//  Recipe from GK10, with full RT
*/
#ifdef RADIATIVE_TRANSFER

struct
{
  double factor;                         /* normalization constant; used to be called fmass */
  double efficiency;                     /* used to be called eps_SFH2 */
  double min_molecular_fraction;         /* used to be called fH2_SFH2 */
  double min_cloud_density;              /* used to be called den_SFH2_eff */
  double max_cloud_density;              /* not used previously; set max_cloud_density = min_cloud_density for recipe 1 of Gnedin et al 2009 */
  double very_high_density;              /* used to be called den_PRIM_eff */
}
sfr = { 0.0, 0.005, 0.1, 50.0, 1.0e99, 1.0e6}; 


void sfr_init()
{
  control_parameter_add3(control_parameter_double,&sfr.efficiency,"sfr:efficiency","sfr.efficiency","sf:recipe=1.efficiency","the efficiency of the star formation law in molecular gas per free-fall time (a-la Krumholz and Tan 2006).");

  control_parameter_add3(control_parameter_double,&sfr.min_molecular_fraction,"sfr:min-molecular-fraction","sfr.min-molecular-fraction","sf:recipe=1.min_molecular_fraction","the minimum molecular (H2) fraction for star formation.");

  control_parameter_add3(control_parameter_double,&sfr.min_cloud_density,"sfr:min-cloud-density","sfr.min-cloud-density","sf:recipe=1:min-cloud-density","the minimum density for computing the free-fall time. The non-zero value of this parameter selects recipes #1 or #2 of Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sfr.max_cloud_density,"sfr:max-cloud-density","sfr.max-cloud-density","sf:recipe=1:max-cloud-density","the maximum density for computing the free-fall time. Setting <sfr:max-cloud-density> = <sfr:min-cloud-density> reduces this recipe to the recipe #1 of Gnedin et al 2009; setting <sfr:max-cloud-density> to a very large number effectively removes this limit and reduces this recipe to the recipe #2 of Gnedin et al 2009; a non-trivial value of <sfr:max-cloud-density> > <sfr:min-cloud-density> makes a recipe not discussed in Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sfr.very_high_density,"sfr:very-high-density","sfr.very-high-density","sf:recipe=1:very-high-density","the minimum density above which all gas is assumed to participate in star formation, irrespectively of its molecular fraction. This can be used to smoothly switch to primordial mode of star formation, when the molecular fraction never exceeds about 0.001.");
}


void sfr_verify()
{
  cart_assert(sfr.efficiency > 0.0);
  cart_assert(sfr.min_molecular_fraction > 0.0);
  cart_assert(!(sfr.min_cloud_density < 0.0));
  cart_assert(!(sfr.max_cloud_density < sfr.min_cloud_density));
  cart_assert(sfr.very_high_density > 0.0);
}


void sfr_setup(int level)
{
  sfr.factor = sfr.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(3*M_PI)); 
}


double sfr_rate(int cell)
{
  double fH2_cell, nH_eff;
  double nH = constants->XH*units->number_density*cell_gas_density(cell);

  fH2_cell = cell_H2_fraction(cell);

  if(nH > sfr.very_high_density) fH2_cell = 1.0;

  nH_eff = max(sfr.min_cloud_density,min(sfr.max_cloud_density,nH));

  if(fH2_cell > sfr.min_molecular_fraction)
    {
      return sfr.factor*fH2_cell*cell_gas_density(cell)*sqrt(nH_eff);
    }
  else
    {
      return 0.0;
    }
}


struct StarFormationRecipe sf_recipe_internal =
{
  "gh10-full",
  sfr_init,
  sfr_verify,
  sfr_setup,
  sfr_rate
};

#else /* RADIATIVE_TRANSFER */

#error "SF Recipe gk10-full only works with RADIATIVE_TRANSFER activated."

#endif /* RADIATIVE_TRANSFER */

#endif /* HYDRO && STARFORM */
