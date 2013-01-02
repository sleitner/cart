#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "starformation_recipe.h"
#include "rand.h"
#include "rt.h"
#include "tree.h"
#include "units.h"


/*
//  Recipe based on Eq. (17) of GK10, with tau_{SF} = const = 1.5Gyr (a-la Genzel et al 2010 & Bigiel et al)
*/
#if defined(RADIATIVE_TRANSFER) && defined(ENRICHMENT)

struct 
{
  double factor;              /* normalization constant; used to be called fmass */
  double efficiency;          /* defined as Sigma_SFR*1.5Gyr/(1.36*Sigma_H2) */
  double variability;
  double min_molecular_fraction;
}
  sfr = { 0.0, 1.0, 0.0, 0.0 };


void sfr_config_init()
{
  control_parameter_add3(control_parameter_double,&sfr.efficiency,"sf:efficiency","sfr.efficiency","sf:recipe=3:efficiency","the relative efficiency of the star formation law (relative to the constant depletion time-scale of 1.5 Gyr).");

  control_parameter_add3(control_parameter_double,&sfr.variability,"sf:variability","sfr.variability","sf:recipe=3:variability","the variability of the efficiency. If <sf:variability> is greater than 0, then it serves as a dispersion of a lognormal distribution with which the actual SF efficiency fluctuates around <sf:efficiency> (i.e. <sf:efficiency> remains the mean of the distribution).");

  control_parameter_add2(control_parameter_double,&sfr.min_molecular_fraction,"sf:min-molecular-fraction","sfr.min-molecular-fraction","the minimum molecular (H2) fraction for star formation.");

}


void sfr_config_verify()
{
  VERIFY(sf:efficiency, sfr.efficiency > 0.0 );
  VERIFY(sf:variability, !(sfr.variability < 0.0) );
  VERIFY(sf:min-molecular-fraction, !(sfr.min_molecular_fraction < 0.0) );
}


void sfr_setup(int level)
{
  sfr.factor = sfr.efficiency*units->time/(1.5*constants->Gyr);
}


double sfr_rate(int cell)
{
#ifdef RT_CHEMISTRY
  /*
  //  Factor 1.36 is what observers add.
  */
  double fracH2 = 2*cell_H2_density(cell)/(0.76*cell_gas_density(cell));
#else 
  double Umw = rtUmw(cell);
  double Dmw = rtDmw2(cell);  /* floor included */
  double Dstar = 8.0e-3*sqrt(0.01+Umw); 
  double g = 0.08/(1+Umw*pow(Dmw/Dstar,6.0));
  double Lambda = log(1+pow(g+pow(Dmw,0.75)*(Umw/15),4.0/7.0));
  double SigmaC = 20*pow(Lambda,4.0/7.0)/(Dmw*sqrt(1+Umw*Dmw*Dmw));

  //double rhoH = cell_HI_density(cell) + 2*cell_H2_density(cell);
  double rhoH = 0.76*cell_gas_density(cell);
  double SigmaH = units->density*units->length/constants->Msun*constants->pc*constants->pc*rhoH*cell_sobolev_length(cell);
  double fac1 = SigmaH/(SigmaH+SigmaC);

  double fracH2 = fac1*fac1;
#endif

  if(fracH2 > sfr.min_molecular_fraction)
    {
      if(sfr.variability > 0.0)
	{
	  return sfr.factor*fracH2*cell_gas_density(cell)*cart_rand_lognormal(sfr.variability);
	}
      else
	{
	  return sfr.factor*fracH2*cell_gas_density(cell);
	}
    }
  else return 0.0;
}


struct StarFormationRecipe sf_recipe_internal =
{
  "linear-in-H2-new",
  sfr_rate,
  sfr_config_init,
  sfr_config_verify,
  sfr_setup
};

#else

#error "SF Recipe line-in-H2 only works with RADIATIVE_TRANSFER and ENRICHMENT activated."

#endif /* RADIATIVE_TRANSFER && ENRICHMENT */

#endif /* HYDRO && STAR_FORMATION */
