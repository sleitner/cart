#include "config.h"
#if defined(HYDRO) && defined(STARFORM)

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "starformation_recipes.h"
#include "rand.h"
#include "rt.h"
#include "tree.h"
#include "units.h"


void sf_recipe0_setup(int level);
void sf_recipe1_setup(int level);
void sf_recipe2_setup(int level);
void sf_recipe3_setup(int level);

double sf_recipe0_rate(int cell);
double sf_recipe1_rate(int cell);
double sf_recipe2_rate(int cell);
double sf_recipe3_rate(int cell);


struct StarFormationRecipe sf_all_recipes[] =
{
  { sf_recipe0_setup, sf_recipe0_rate, 0, "hart" },
  { sf_recipe1_setup, sf_recipe1_rate, 1, "gk10-full", },
  { sf_recipe2_setup, sf_recipe2_rate, 2, "gk10-lite", },
  { sf_recipe3_setup, sf_recipe3_rate, 3, "linear" }
};


const struct StarFormationRecipe *sf_recipe = sf_all_recipes + 0; //NULL;  /* there is no default recipe */


double sf_factor = 0.0;  /* normalization constant; used to be called fmass */


/*
//  Recipe #0
//  --------------------------------------------------------
*/
struct 
{
  double slope;               /* used to be called alpha_SF */
  double efficiency;          /* used to be called eps_SF */
}
sf_recipe0 = { 1.5, 1.5 };

void sf_recipe0_setup(int level)
{
  sf_factor = sf_recipe0.efficiency*units->time/(4.0e9*constants->yr)*pow(units->density*pow(constants->Mpc,3.0)/(1.0e16*constants->Msun),sf_recipe0.slope-1);
}

double sf_recipe0_rate(int cell)
{
  return sf_factor*pow(cell_gas_density(cell),sf_recipe0.slope);
}


/*
//  Recipe #1
//  --------------------------------------------------------
*/
struct
{
  double efficiency;                     /* used to be called eps_SFH2 */
  double min_molecular_fraction;         /* used to be called fH2_SFH2 */
  double min_cloud_density;              /* used to be called den_SFH2_eff */
  double max_cloud_density;              /* not used previously; set max_cloud_density = min_cloud_density for recipe 1 of Gnedin et al 2009 */
  double very_high_density;              /* used to be called den_PRIM_eff */
}
sf_recipe1 = { 0.005, 0.1, 50.0, 1.0e99, 1.0e6}; 

void sf_recipe1_setup(int level)
{
  sf_factor = sf_recipe1.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(3*M_PI)); 
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


/*
//  Recipe #2
//  --------------------------------------------------------
*/
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

void sf_recipe2_setup(int level)
{
  sf_factor = sf_recipe2.efficiency*units->time*sqrt(32*constants->G*constants->XH*constants->mp/(3*M_PI)); 
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
  cart_error("SF Recipe #2: only works with ENRICH activated.");
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

/*
//  Recipe #3
//  --------------------------------------------------------
*/
struct 
{
  double efficiency;          /* defined as Sigma_SFR*1.5Gyr/(1.36*Sigma_H2) */
  double variability;
}
sf_recipe3 = { 1.0, 0.0 };

void sf_recipe3_setup(int level)
{
  sf_factor = sf_recipe3.efficiency*units->time/(1.5*constants->Gyr);
}

double sf_recipe3_rate(int cell)
{
#if defined(RADIATIVE_TRANSFER) && defined(ENRICH)
#ifdef RT_CHEMISTRY
  /*
  //  Factor 1.36 is what observers add.
  */
  double rhoH2 = 1.36*2*cell_H2_density(cell);
#else 
  double Umw = rtUmw(cell);
  double Dmw = rtDmw2(cell);  /* floor included */
  double alpha = 5*(Umw/2)/(1+pow(Umw/2,2.0)); 
  double Dstar = 1.5e-3*log(1+pow(3*Umw,1.7)); 
  double s = 0.04/(Dstar+Dmw); 
  double g = (1+alpha*s+s*s)/(1+s); 
  double Lambda = log(1+g*pow(Dmw,3.0/7.0)*pow(Umw/15,4.0/7.0));
  double SigmaC = 20*pow(Lambda,4.0/7.0)/(Dmw*sqrt(1+Umw*Dmw*Dmw));

  double rhoH = cell_HI_density(cell) + 2*cell_H2_density(cell);
  double SigmaH = units->density*units->length/constants->Msun*constants->pc*constants->pc*rhoH*cell_sobolev_length(cell);
  double fac1 = SigmaH/(SigmaH+SigmaC);

  double rhoH2 = 1.36*rhoH*fac1*fac1;
#endif

  if(sf_recipe3.variability > 1.0)
    {
      if(sf_recipe3.variability*cart_rand() < 1.0)
	{
	  return sf_recipe3.variability*sf_factor*rhoH2;
	}
      else return 0.0;
    }
  else
    {
      return sf_factor*rhoH2;
    }

#else
  cart_error("SF Recipe #3: only works with RADIATIVE_TRANSFER and ENRICH activated.");
  return 0.0;
#endif /* RADIATIVE_TRANSFER && ENRICH */
}


/*
//  Configuration
*/
void control_parameter_set_recipe(const char *value, void *ptr, int ind)
{
  const int n = sizeof(sf_all_recipes)/sizeof(struct StarFormationRecipe);
  int i, id;
  char c;

  sf_recipe = NULL;

  /*
  //  First, try to read as an int
  */
  if(sscanf(value,"%d%c",&id,&c) == 1)
    {
      for(i=0; i<n; i++)
	{
	  if(sf_all_recipes[i].id == id)
	    {
	      sf_recipe = sf_all_recipes + i;
	      break;
	    }
	}
    }
  else
    {
      /*
      //  Ok, int read failed, try a string
      */
      for(i=0; i<n; i++)
	{
	  if(strstr(sf_all_recipes[i].name,value) == sf_all_recipes[i].name)
	    {
	      if(sf_recipe == NULL)
		{
		  sf_recipe = sf_all_recipes + i;
		}
	      else
		{
		  cart_error("String <%s> is an ambigous name for a SF recipe.",value);
		}
	    }
        }
    }

  if(sf_recipe == NULL)
    {
      cart_debug("String '%s' is not a valid specifier for a SF recipe. Valid ids are:",value);
      for(i=0; i<n; i++) cart_debug("%s or %d",sf_all_recipes[i].name,sf_all_recipes[i].id);
      cart_error("ART is terminating.");
    }
}


void control_parameter_list_recipe(FILE *stream, const void *ptr)
{
  cart_assert(sf_recipe != NULL);

  fprintf(stream,"%s (or id=%d)",sf_recipe->name,sf_recipe->id);
}


void config_init_star_formation_recipes()
{
  ControlParameterOps control_parameter_recipe = { control_parameter_set_recipe, control_parameter_list_recipe };

  control_parameter_add2(control_parameter_recipe,&sf_recipe,"sf:recipe","sf_recipe","recipe for star formation. Available recipes: \n   'hart' or '0' (oldstyle HART recipe),\n   'gk10-full' or '1' (Gnedin et al 2009 recipes),\n   'gk10-lite' or '2' (Gnedin & Kravtsov 2010 eq 6 recipe),\n   'linear' or '3' (linear in H2 density, without H2 falls back onto Gnedin & Kravtsov 2010 eq 17 recipe).");

  /*
  //  Recipe 0 (oldstyle HART recipe )
  */
  control_parameter_add4(control_parameter_double,&sf_recipe0.slope,"sf:recipe=<hart>:slope","sf:recipe=0:slope","sf_recipe0.slope","alpha_sf","the slope of the star formation law with gas density (see HART documentation for exact definition).");

  control_parameter_add4(control_parameter_double,&sf_recipe0.efficiency,"sf:recipe=<hart>:efficiency","sf:recipe=0:efficiency","sf_recipe0.efficiency","eps_sf","the relative efficiency of the star formation law (see HART documentation for exact definition).");

  /*
  //  Recipe 1 (combines recipes 1-3 of Gnedin et al 2009)
  */
  control_parameter_add3(control_parameter_double,&sf_recipe1.efficiency,"sf:recipe=<gk10-full>:efficiency","sf:recipe=1:efficiency","sf_recipe1.efficiency","the efficiency of the star formation law in molecular gas per free-fall time (a-la Krumholz and Tan 2006).");

  control_parameter_add3(control_parameter_double,&sf_recipe1.min_molecular_fraction,"sf:recipe=<gk10-full>:min-molecular-fraction","sf:recipe=1:min-molecular-fraction","sf_recipe1.min_molecular_fraction","the minimum molecular (H2) fraction for star formation.");

  control_parameter_add3(control_parameter_double,&sf_recipe1.min_cloud_density,"sf:recipe=<gk10-full>:min-cloud-density","sf:recipe=1:min-cloud-density","sf_recipe1.min_cloud_density","the minimum density for computing the free-fall time. The non-zero value of this parameter selects recipes #1 or #2 of Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sf_recipe1.max_cloud_density,"sf:recipe=<gk10-full>:max-cloud-density","sf:recipe=1:max-cloud-density","sf_recipe1.max_cloud_density","the maximum density for computing the free-fall time. Setting <sf:recipe=1:max-cloud-density> = <sf:recipe=1:min-cloud-density> reduces to the recipe #1 of Gnedin et al 2009; setting <sf:recipe=1:max-cloud-density> to a very large number effectively removes this limit and reduces to the recipe #2 of Gnedin et al 2009; a non-trivial value of <sf:recipe=1:max-cloud-density> > <sf:recipe=1:min-cloud-density> makes a recipe now discussed in Gnedin et al 2009.");

  control_parameter_add3(control_parameter_double,&sf_recipe1.very_high_density,"sf:recipe=<gk10-full>:very-high-density","sf:recipe=1:very-high-density","sf_recipe1.very_high_density","the minimum density above which all gas is assumed to participate in star formation, irrespectively of its molecular fraction. This can be used to smoothly switch to primordial mode of star formation, when the molecular fraction never exceeds about 0.001.");

  /*
  //  Recipe 2 (Equation 6 of Gnedin & Kravtsov 2010)
  */
  control_parameter_add3(control_parameter_double,&sf_recipe2.efficiency,"sf:recipe=<gk10-lite>:efficiency","sf:recipe=2:efficiency","sf_recipe2.efficiency","the efficiency of the star formation law in molecular gas per free-fall time (a-la Krumholz and Tan 2006).");

  control_parameter_add3(control_parameter_double,&sf_recipe2.min_molecular_fraction,"sf:recipe=<gk10-lite>:min-molecular-fraction","sf:recipe=2:min-molecular-fraction","sf_recipe2.min_molecular_fraction","the minimum molecular (H2) fraction for star formation.");

  control_parameter_add3(control_parameter_double,&sf_recipe2.min_cloud_density,"sf:recipe=<gk10-lite>:min-cloud-density","sf:recipe=2:min-cloud-density","sf_recipe2.min_cloud_density","the minimum density for computing the free-fall time. ");

  control_parameter_add3(control_parameter_double,&sf_recipe2.max_cloud_density,"sf:recipe=<gk10-lite>:max-cloud-density","sf:recipe=2:max-cloud-density","sf_recipe2.max_cloud_density","the maximum density for computing the free-fall time.");

  control_parameter_add3(control_parameter_double,&sf_recipe2.very_high_density,"sf:recipe=<gk10-lite>:very-high-density","sf:recipe=2:very-high-density","sf_recipe2.very_high_density","the minimum density above which all gas is assumed to participate in star formation, irrespectively of its molecular fraction. This can be used to smoothly switch to primordial mode of star formation, when the molecular fraction never exceeds about 0.001.");

  control_parameter_add3(control_parameter_double,&sf_recipe2.U_MW,"sf:recipe=<gk10-lite>:U_MW","sf:recipe=2:U_MW","sf_recipe2.U_MW","If positive, choose a constant value for the U_MW in Gnedin & Kravtsov 2010 eq 6. If negative you are selecting an estimate for U_MW based on the SFR smoothed on some scale... under development");

  control_parameter_add3(control_parameter_double,&sf_recipe2.D_MW,"sf:recipe=<gk10-lite>:D_MW","sf:recipe=2:D_MW","sf_recipe2.D_MW","If positive, choose a constant value for the D_MW in Gnedin & Kravtsov 2010 eq 6. If negative or zero you are selecting an estimate for D_MW based on the cell metallicity.");

  /*
  //  Recipe 3 (Eq. (17) GK10)
  */
  control_parameter_add3(control_parameter_double,&sf_recipe3.efficiency,"sf:recipe=<linear>:efficiency","sf:recipe=3:efficiency","sf_recipe3.efficiency","the relative efficiency of the star formation law (see GF10 documentation for exact definition).");

  control_parameter_add3(control_parameter_double,&sf_recipe3.variability,"sf:recipe=<linear>:variability","sf:recipe=3:variability","sf_recipe3.variability","the variability of the efficiency.");

}


void config_verify_star_formation_recipes()
{
  cart_assert(sf_recipe != NULL);

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

  /*
  //  Recipe 2 (Gnedin & Kravtsov 2010 eq 17)
  */
  cart_assert(sf_recipe3.efficiency > 0.0);
  cart_assert(sf_recipe3.variability==0.0 || sf_recipe3.variability>=1.0);
}

#endif /* HYDRO && STARFORM */
