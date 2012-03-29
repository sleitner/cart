#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "imf.h"


typedef double(*fimf)(double);

double f_IMF_Salpeter(double m);
double f_IMF_MillerScalo(double m);
double f_IMF_Chabrier(double m);
double f_IMF_Kroupa(double m);


struct IMF_t
{
  char* name;
  fimf  f;
}
const IMF_fun[] = { 
  { "Salpeter", f_IMF_Salpeter }, 
  { "Miller-Scalo", f_IMF_MillerScalo }, 
  { "Chabrier", f_IMF_Chabrier }, 
  { "Kroupa", f_IMF_Kroupa } 
};

const int num_imfs = sizeof(IMF_fun)/sizeof(struct IMF_t);


struct InitialMassFunction imf_internal = { 1, 0.1, 100.0, 8.0, 3.0, 8.0 };
const struct InitialMassFunction *imf = &imf_internal;


void control_parameter_set_imf(const char *value, void *ptr, int ind)
{
  int i;

  imf_internal.type = -1;

  for(i=0; i<num_imfs; i++)
    {
      if(strcmp(value,IMF_fun[i].name) == 0)
	{
	  imf_internal.type = i;
	}
    }

  if(imf_internal.type < 0)
    {
      cart_debug("String '%s' is not a valid name of IMF. Valid names are:",value);
      for(i=0; i<num_imfs; i++) cart_debug("%s",IMF_fun[i].name);
      cart_error("ART is terminating.");
    }
}


void control_parameter_list_imf(FILE *stream, const void *ptr)
{
  control_parameter_list_string(stream,IMF_fun[imf_internal.type].name);
}


void config_init_imf()
{
  ControlParameterOps control_parameter_imf = { control_parameter_set_imf, control_parameter_list_imf };

  control_parameter_add2(control_parameter_imf,&imf_internal.type,"imf","imf.type","the name of the IMF. Valid names are 'Salpeter', 'Miller-Scalo', and 'Chabrier'.");

  control_parameter_add3(control_parameter_double,&imf_internal.min_mass,"imf:min-mass","imf.min_mass","am_stl","the minimum stellar mass in the IMF model.");

  control_parameter_add3(control_parameter_double,&imf_internal.max_mass,"imf:max-mass","imf.max_mass","am_stu","the maximum stellar mass in the IMF model.");

  control_parameter_add3(control_parameter_double,&imf_internal.min_SNII_mass,"imf:min-SNII-mass","imf.min_SNII_mass","am_snii","the minimum mass of stars that explode as type II supernovae.");

  control_parameter_add3(control_parameter_double,&imf_internal.min_SNIa_mass,"imf:min-SNIa-mass","imf.min_SNIa_mass","am_snia1","the minimum mass of stars that explode as type Ia supernovae.");

  control_parameter_add3(control_parameter_double,&imf_internal.max_SNIa_mass,"imf:max-SNIa-mass","imf.max_SNIa_mass","am_snia2","the maximum mass of stars that explode as type Ia supernovae.");

}


void config_verify_imf()
{
  /*
  //  IMF
  */
  cart_assert(imf_internal.type>=0 && imf_internal.type<num_imfs);

  cart_assert(imf_internal.min_mass > 0.0);

  cart_assert(imf_internal.max_mass > 0.0);

  cart_assert(imf_internal.min_SNII_mass > 1.0);

  cart_assert(imf_internal.min_SNIa_mass > 1.0);

  cart_assert(imf_internal.max_SNIa_mass > imf_internal.min_SNIa_mass);

}




double f_IMF_Salpeter( double m )
{
  return pow( m, -2.35 );
}


double f_IMF_MillerScalo( double m )
{
  /* Miller-Scalo (1979, ApJS 41, 513, eq. 30, Table 7) IMF */
  const double C_1 =  1.09;
  const double C_2 = -1.02;
  return exp( -C_1 * ( (log10(m) - C_2)*(log10(m) - C_2) ) ) / m;
}


double f_IMF_Chabrier( double m )
{
  /* Chabrier, G. (2001, ApJ 554, 1274) */
  const double m0_Ch = 716.4;
  const double beta_Ch = 0.25;
  const double alpha_Ch = -3.3;
  return exp( -pow( m0_Ch/m, beta_Ch) ) * pow(m,alpha_Ch);
}


double f_IMF_Kroupa( double m )
{
  /* Kroupa, P. (2007astro.ph..3124K) */
  if(m < 0.5)
    return pow(m/0.08,-1.3);
  else
    return pow(0.5/0.08,-1.3)*pow(m/0.5,-2.3);
}
  

double f_IMF( double mstar )
{
  return IMF_fun[imf_internal.type].f(mstar);
}


double fm_IMF( double mstar )
{
  return mstar * f_IMF(mstar);
}


double fmz_IMF( double mstar )
{
  return mstar * f_IMF(mstar) * min( 0.2, max( 0.01*mstar - 0.06, 1e-20 ) );
}

#endif /* STAR_FORMATION */
