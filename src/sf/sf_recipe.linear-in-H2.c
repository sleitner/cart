#include "config.h"
#if defined(HYDRO) && defined(STARFORM)

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
#if defined(RADIATIVE_TRANSFER) && defined(ENRICH)

struct 
{
  double factor;              /* normalization constant; used to be called fmass */
  double efficiency;          /* defined as Sigma_SFR*1.5Gyr/(1.36*Sigma_H2) */
  double variability;
}
sfr = { 0.0, 1.0, 0.0 };


void sfr_init()
{
  control_parameter_add3(control_parameter_double,&sfr.efficiency,"sfr:efficiency","sfr.efficiency","sf:recipe=3:efficiency","the relative efficiency of the star formation law (relative to the constant depletion time-scale of 1.5 Gyr).");

  control_parameter_add3(control_parameter_double,&sfr.variability,"sfr:variability","sfr.variability","sf:recipe=3:variability","the variability of the efficiency. If <sfr:variability> is greater than 0, then it serves as a dispersion of a lognormal distribution with which the actual SF efficiency fluctuates around <sfr:efficiency> (i.e. <sfr:efficiency> remains the mean of the distribution).");
}


void sfr_verify()
{
  cart_assert(sfr.efficiency > 0.0);
  cart_assert(!(sfr.variability < 0.0));
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

  //double rhoH = cell_HI_density(cell) + 2*cell_H2_density(cell);
  double rhoH = 0.76*cell_gas_density(cell);
  double SigmaH = units->density*units->length/constants->Msun*constants->pc*constants->pc*rhoH*cell_sobolev_length(cell);
  double fac1 = SigmaH/(SigmaH+SigmaC);

  double rhoH2 = 1.36*rhoH*fac1*fac1;
#endif

  if(sfr.variability > 0.0)
    {
      return sfr.factor*rhoH2*cart_rand_lognormal(sfr.variability);
    }
  else
    {
      return sfr.factor*rhoH2;
    }
}


struct StarFormationRecipe sf_recipe_internal =
{
  "linear-in-H2",
  sfr_init,
  sfr_verify,
  sfr_setup,
  sfr_rate
};

#else

#error "SF Recipe line-in-H2 only works with RADIATIVE_TRANSFER and ENRICH activated."

#endif /* RADIATIVE_TRANSFER && ENRICH */

#endif /* HYDRO && STARFORM */
