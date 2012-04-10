#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_recipe.h"
#include "tree.h"
#include "units.h"


/*
//  Recipe from GK10, without full RT (using equation 6)
*/
struct
{
  double factor;                         /* normalization constant; used to be called fmass */
  double efficiency;                     /* used to be called eps_SFH2 */
  double min_molecular_fraction;         /* used to be called fH2_SFH2 */
  double min_cloud_density;              /* used to be called den_SFH2_eff */
  double max_cloud_density;              /* not used previously; set max_cloud_density = min_cloud_density for recipe 1 of Gnedin et al 2009 */
  double very_high_density;              /* used to be called den_PRIM_eff */
  double U_MW;                           /* not used previously; if <=0: set a fit to UV, if positive: set it as a constant */
  double D_MW;                           /* not used previously; if <=0: set D_MW to metallicity, if positive: set it as a constant */
}
sfr = { 0.0, 0.005, 0.1, 50.0, 1.0e99, 1.0e6, 1.0, -1 }; 


void sfr_config_init()
{
  control_parameter_add3(control_parameter_double,&sfr.efficiency,"sf:efficiency","sfr.efficiency","sf:recipe=1.efficiency","the efficiency of the star formation law in molecular gas per free-fall time (a-la Krumholz and Tan 2006).");

  control_parameter_add3(control_parameter_double,&sfr.min_molecular_fraction,"sf:min-molecular-fraction","sfr.min-molecular-fraction","sf:recipe=1.min_molecular_fraction","the minimum molecular (H2) fraction for star formation.");

  control_parameter_add3(control_parameter_double,&sfr.min_cloud_density,"sf:min-cloud-density","sfr.min-cloud-density","sf:recipe=1:min-cloud-density","the minimum density for computing the free-fall time. The non-zero value of this parameter selects recipes #1 or #2 of Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sfr.max_cloud_density,"sf:max-cloud-density","sfr.max-cloud-density","sf:recipe=1:max-cloud-density","the maximum density for computing the free-fall time. Setting <sf:max-cloud-density> = <sf:min-cloud-density> reduces this recipe to the recipe #1 of Gnedin et al 2009; setting <sf:max-cloud-density> to a very large number effectively removes this limit and reduces this recipe to the recipe #2 of Gnedin et al 2009; a non-trivial value of <sf:max-cloud-density> > <sf:min-cloud-density> makes a recipe not discussed in Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sfr.very_high_density,"sf:very-high-density","sfr.very-high-density","sf:recipe=1:very-high-density","the minimum density above which all gas is assumed to participate in star formation, irrespectively of its molecular fraction. This can be used to smoothly switch to primordial mode of star formation, when the molecular fraction never exceeds about 0.001.");

  control_parameter_add3(control_parameter_double,&sfr.U_MW,"sf:U_MW","sfr.U_MW","sf:recipe=2:U_MW","If positive, choose a constant value for the U_MW in Gnedin & Kravtsov 2010 eq 6. If -1, you are selecting an estimate for U_MW based on the SFR smoothed on some scale... (under development).");

  control_parameter_add3(control_parameter_double,&sfr.D_MW,"sf:D_MW","sfr.D_MW","sf:recipe=2:D_MW","If positive, choose a constant value for the D_MW in Gnedin & Kravtsov 2010 eq 6. If -1, you are selecting an estimate for D_MW based on the cell metallicity.");
}


void sfr_config_verify()
{
  VERIFY(sf:efficiency, sfr.efficiency > 0.0 );
  VERIFY(sf:min-molecular-fraction, sfr.min_molecular_fraction > 0.0 );
  VERIFY(sf:min-cloud-density, !(sfr.min_cloud_density < 0.0) );
  VERIFY(sf:max-cloud-density, !(sfr.max_cloud_density < sfr.min_cloud_density) );
  VERIFY(sf:very-high-density, sfr.very_high_density > 0.0 );
  VERIFY(sf:U_MW, (sfr.U_MW==-1 || sfr.U_MW>0.0) );
  VERIFY(sf:D_MW, (sfr.D_MW==-1 || sfr.D_MW>0.0) );
}


void sfr_setup(int level)
{
  sfr.factor = sfr.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(3*M_PI)); 
}


double sfr_rate(int cell)
{
  double fH2_cell, nH_eff;
  double nH = constants->XH*units->number_density*cell_gas_density(cell);
  
  double zSol_cell;
  double D_MW,U_MW; 
  double Dstar,g,s,x,alpha,lambda;
  double nstar = 25; 

  if(sfr.U_MW > 0)
    {
      U_MW = sfr.U_MW;
    }
  else
    {
      /* UV <= unknown  esp at high z */
      cart_error("SF Recipe gk10-lite: fit U_MW to some local SF estimate. This has not been done yet.");
    }

  if(sfr.D_MW > 0)
    {
      D_MW = sfr.D_MW;
    }
  else
    {
#ifdef ENRICHMENT
      zSol_cell = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
      D_MW = max(1.0e-3,zSol_cell);
#else
      cart_error("SF Recipe gk10-lite with D_MW<0 only works with ENRICHMENT activated.");
      D_MW = 1.0e-3;
#endif /* ENRICHMENT */
    }
  
  alpha = 5 * ( U_MW/2. )/( 1 + pow(U_MW/2.,2)  ); 
  Dstar = 1.5e-3*log( 1 + pow(3*U_MW,1.7) ); 
  s = 0.04 /( Dstar + D_MW ); 
  g = (1 + alpha*s + s*s )/( 1 + s ); 
  lambda = log( 1 + g*pow(D_MW,3./7.)*pow(U_MW/15.,4./7.)  ); 
  x = pow(lambda,3./7.) * log( D_MW*nH / (lambda*nstar) ) ; 
  
  fH2_cell = 1 / ( 1 + exp(-4*x - 3*pow(x,3)) ); 
  if( fH2_cell < 0.1 ){
    x=x/pow(g,0.25);
    fH2_cell = 1 / ( 1 + exp(-4*x - 3*pow(x,3)) );
  }
  
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
  "gk10-lite",
  sfr_rate,
  sfr_config_init,
  sfr_config_verify,
  sfr_setup
};

#endif /* HYDRO && STAR_FORMATION */
