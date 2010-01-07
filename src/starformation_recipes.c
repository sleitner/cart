#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "config.h"
#include "control_parameter.h"
#include "starformation.h"
#include "starformation_recipes.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"


extern int sf_recipe;


StarFormationRate sf_rate = NULL;
double sf_factor = 0.0;  /* normalization constant; used to be called fmass */

struct 
{
  double slope;               /* used to be called alpha_SF */
  double efficiency;          /* used to be called eps_SF */
}
sf_recipe0 = { 1.5, 1.5 };

struct
{
  double efficiency;                     /* used to be called eps_SFH2 */
  double min_molecular_fraction;         /* used to be called fH2_SFH2 */
  double min_cloud_density;              /* used to be called den_SFH2_eff */
  double max_cloud_density;              /* not used previously; set max_cloud_density = min_cloud_density for recipe 1 of Gnedin et al 2009 */
  double very_high_density;              /* used to be called den_PRIM_eff */
}
sf_recipe1 = { 0.005, 0.1, 50.0, 1.0e99, 1.0e6 };


void config_init_star_formation_recipes()
{
  /*
  //  Recipe 0 (oldstyle HART recipe )
  */
  control_parameter_add3(control_parameter_double,&sf_recipe0.slope,"sf:recipe=0:slope","sf_recipe0.slope","alpha_sf","the slope of the star formation law with gas density (see HART documentation for exact definition).");

  control_parameter_add3(control_parameter_double,&sf_recipe0.efficiency,"sf:recipe=0:efficiency","sf_recipe0.efficiency","eps_sf","the relative efficiency of the star formation law (see HART documentation for exact definition).");

  /*
  //  Recipe 1 (combines recipes 1-3 of Gnedin et al 2009)
  */
  control_parameter_add3(control_parameter_double,&sf_recipe1.efficiency,"sf:recipe=1:efficiency","sf_recipe1.efficiency","eps_sfh2","the efficiency of the star formation law in molecular gas per free-fall time (a-la Krumholz and Tan 2006).");

  control_parameter_add3(control_parameter_double,&sf_recipe1.min_molecular_fraction,"sf:recipe=1:min-molecular-fraction","sf_recipe1.min_molecular_fraction","fh2_sfh2","the minimum molecular (H2) fraction for star formation.");

  control_parameter_add3(control_parameter_double,&sf_recipe1.min_cloud_density,"sf:recipe=1:min-cloud-density","sf_recipe1.min_cloud_density","den_sfh2_eff","the minimum density for computing the free-fall time. The non-zero value of this parameter selects recipes #1 or #2 of Gnedin et al 2009.");

  control_parameter_add2(control_parameter_double,&sf_recipe1.max_cloud_density,"sf:recipe=1:max-cloud-density","sf_recipe1.max_cloud_density","the maximum density for computing the free-fall time. Setting <sf:recipe=1:max-cloud-density> = <sf:recipe=1:min-cloud-density> reduces to the recipe #1 of Gnedin et al 2009; setting <sf:recipe=1:max-cloud-density> to a very large number effectively removes this limit and reduces to the recipe #2 of Gnedin et al 2009; a non-trivial value of <sf:recipe=1:max-cloud-density> > <sf:recipe=1:min-cloud-density> makes a recipe now discussed in Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sf_recipe1.very_high_density,"sf:recipe=1:very-high-density","sf_recipe1.very_high_density","den_prim_eff","the minimum density above which all gas is assumed to participate in star formation, irrespectively of its molecular fraction. This can be used to smoothly switch to primordial mode of star formation, when the molecular fraction never exceeds about 0.001.");
}


void config_verify_star_formation_recipes()
{
  /*
  //  Recipe 0 (oldstyle HART recipe )
  */
  cart_assert(sf_recipe0.slope > 0.0);

  cart_assert(sf_recipe0.efficiency > 0.0);

  /*
  //  Recipe 1 (combines recipes 1-3 of Gnedin et al 2009)
  */
  cart_assert(sf_recipe1.efficiency > 0.0);

  cart_assert(sf_recipe1.min_molecular_fraction > 0.0);

  cart_assert(!(sf_recipe1.min_cloud_density < 0.0));

  cart_assert(!(sf_recipe1.max_cloud_density < sf_recipe1.min_cloud_density));

  cart_assert(sf_recipe1.very_high_density > 0.0);
}


#ifdef HYDRO

double sf_recipe0_rate(int cell)
{
  return sf_factor*pow(cell_gas_density(cell),sf_recipe0.slope);
}


double sf_recipe1_rate(int cell)
{
  double fH2_cell, nH_eff;
  double nH = constants->XH*units->number_density*cell_gas_density(cell);

#ifdef RADIATIVE_TRANSFER

  fH2_cell = 2*cell_H2_density(cell)/(2*cell_H2_density(cell)+cell_HI_density(cell));

#else /* RADIATIVE_TRANSFER */
  double zSol_cell;

#ifdef ENRICH
  zSol_cell = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  zSol_cell = 0.0;
#endif /* ENRICH */

  fH2_cell = (max(1.0e-3,zSol_cell)*nH > 30.0) ? 1.0 : 0.0;

#endif /* RADIATIVE_TRANSFER */

  if(nH > sf_recipe1.very_high_density) fH2_cell = 1.0;

  nH_eff = max(sf_recipe1.min_cloud_density,min(sf_recipe1.max_cloud_density,nH));

  if(fH2_cell > sf_recipe1.min_molecular_fraction)
    {
      return sf_factor*fH2_cell*cell_gas_density(cell)*sqrt(nH_eff);
    }
  else
    {
      return 0.0;
    }
}


void setup_star_formation_recipes(int level)
{
  switch(sf_recipe)
    {
    case 0:
      {
	sf_factor = sf_recipe0.efficiency*units->time/(4.0e9*constants->yr)*pow(units->density*pow(constants->Mpc,3.0)/(1.0e16*constants->Msun),sf_recipe0.slope-1);
	sf_rate = sf_recipe0_rate;
	break;
      }
    case 1:
      {
	sf_factor = sf_recipe1.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(2*M_PI));
	sf_rate = sf_recipe1_rate;
	break;
      }
    default:
      {
	sf_rate = NULL;
	cart_error("Invalid sf_recipe value: %d",sf_recipe);
      }
    }
}


#endif /* HYDRO */
#endif /* STARFORM */

