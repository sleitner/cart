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
sf_recipe1 = { 0.005, 0.1, 50.0, 1.0e99, 1.0e6}; 

struct
{
  double efficiency;                     /* used to be called eps_SFH2 */
  double min_molecular_fraction;         /* used to be called fH2_SFH2 */
  double min_cloud_density;              /* used to be called den_SFH2_eff */
  double max_cloud_density;              /* not used previously;  */
  double very_high_density;              /* used to be called den_PRIM_eff */
  double U_MW;                           /* not used previously; if <=0: set a fit to UV, if positive: set it as a constant */
  double D_MW;                           /* not used previously; if <=0: set D_MW to metallicity, if positive: set it as a constant */
}
sf_recipe2 = { 0.005, 0.1, 50.0, 1.0e99, 1.0e6, 1.0, -1}; 

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

  /*
  //  Recipe 2 (Equation 6 of Gnedin & Kravtsov 2010)
  */
  control_parameter_add3(control_parameter_double,&sf_recipe2.efficiency,"sf:recipe=2:efficiency","sf_recipe2.efficiency","eps_sfh2","the efficiency of the star formation law in molecular gas per free-fall time (a-la Krumholz and Tan 2006).");

  control_parameter_add3(control_parameter_double,&sf_recipe2.min_molecular_fraction,"sf:recipe=2:min-molecular-fraction","sf_recipe2.min_molecular_fraction","fh2_sfh2","the minimum molecular (H2) fraction for star formation.");

  control_parameter_add3(control_parameter_double,&sf_recipe2.min_cloud_density,"sf:recipe=2:min-cloud-density","sf_recipe2.min_cloud_density","den_sfh2_eff","the minimum density for computing the free-fall time. ");

  control_parameter_add2(control_parameter_double,&sf_recipe2.max_cloud_density,"sf:recipe=2:max-cloud-density","sf_recipe2.max_cloud_density","the maximum density for computing the free-fall time.");

  control_parameter_add3(control_parameter_double,&sf_recipe2.very_high_density,"sf:recipe=2:very-high-density","sf_recipe2.very_high_density","den_prim_eff","the minimum density above which all gas is assumed to participate in star formation, irrespectively of its molecular fraction. This can be used to smoothly switch to primordial mode of star formation, when the molecular fraction never exceeds about 0.001.");

  control_parameter_add2(control_parameter_double,&sf_recipe2.U_MW,"sf:recipe=2:U_MW","sf_recipe2.U_MW","If positive, choose a constant value for the U_MW in Gnedin & Kravtsov 2010 eq 6. If negative you are selecting an estimate for U_MW based on the SFR smoothed on some scale... under development");

  control_parameter_add2(control_parameter_double,&sf_recipe2.D_MW,"sf:recipe=2:D_MW","sf_recipe2.D_MW","If positive, choose a constant value for the D_MW in Gnedin & Kravtsov 2010 eq 6. If negative or zero you are selecting an estimate for D_MW based on the cell metallicity.");
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

  /*
  //  Recipe 2 (Gnedin & Kravtsov 2010 eq 6)
  */
  cart_assert(sf_recipe2.efficiency > 0.0);

  cart_assert(sf_recipe2.min_molecular_fraction > 0.0);

  cart_assert(!(sf_recipe2.min_cloud_density < 0.0));

  cart_assert(!(sf_recipe2.max_cloud_density < sf_recipe2.min_cloud_density));

  cart_assert(sf_recipe2.very_high_density > 0.0);
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

  fH2_cell = cell_H2_fraction(cell);

#else /* RADIATIVE_TRANSFER */

  cart_error("SF Recipe #1: only works with RADIATIVE_TRANSFER activated.");

  fH2_cell = 0.0;

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


double sf_recipe2_rate(int cell)
{
  double fH2_cell, nH_eff;
  double nH = constants->XH*units->number_density*cell_gas_density(cell);
  
  double zSol_cell;
  double D_MW,U_MW; 
  double Dstar,g,s,x,alpha,lambda;
  double nstar = 25; 

#ifdef ENRICH
  zSol_cell = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  cart_error("ERROR: Need enrichment for SF Recipe #2");
  zSol_cell = 0.0;
#endif /* ENRICH */

  zSol_cell = max(1.0e-3,zSol_cell);
  if( sf_recipe2.U_MW >0 ){
    U_MW = sf_recipe2.U_MW ;
  }else{
    //UV <= unknown  esp at high z
    cart_error("SF Recipe #2: fit U_MW to some local SF estimate. This has not been done yet.");
  }
  if( sf_recipe2.D_MW >0 ){
    D_MW = sf_recipe2.D_MW ;
  }else{
    D_MW = zSol_cell ; 
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
  
  if(nH > sf_recipe2.very_high_density) fH2_cell = 1.0;
  
  nH_eff = max(sf_recipe2.min_cloud_density,min(sf_recipe2.max_cloud_density,nH));
  
  if(fH2_cell > sf_recipe2.min_molecular_fraction)
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
	sf_factor = sf_recipe1.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(3*M_PI)); 
	sf_rate = sf_recipe1_rate;
	break;
      }
    case 2:
      {
	sf_factor = sf_recipe2.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(3*M_PI)); 
	sf_rate = sf_recipe2_rate;
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

