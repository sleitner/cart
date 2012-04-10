#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "imf.h"


double imf_f_Salpeter(double m);
double imf_f_MillerScalo(double m);
double imf_f_Chabrier(double m);
double imf_f_Kroupa(double m);

double imf_fm( double mstar );
double imf_fmz( double mstar );


struct imf_t
{
  const char* name;
  double (*f)(double m);
}
const imf_list[] = { 
  { "Salpeter", imf_f_Salpeter }, 
  { "Miller-Scalo", imf_f_MillerScalo }, 
  { "Chabrier", imf_f_Chabrier }, 
  { "Kroupa", imf_f_Kroupa } 
};

const int num_imfs = sizeof(imf_list)/sizeof(struct imf_t);


struct InitialMassFunction imf_internal = { NULL, 0.1, 100.0, NULL, imf_fm, imf_fmz };
const struct InitialMassFunction *imf = &imf_internal;


void control_parameter_set_imf(const char *value, void *ptr, int ind)
{
  int i;

  imf_internal.name = NULL;

  for(i=0; i<num_imfs; i++)
    {
      if(strcmp(value,imf_list[i].name) == 0)
	{
	  imf_internal.name = imf_list[i].name;
	  imf_internal.f = imf_list[i].f;
	}
    }

  if(imf_internal.name == NULL)
    {
      cart_debug("String '%s' is not a valid name of IMF. Valid names are:",value);
      for(i=0; i<num_imfs; i++) cart_debug("%s",imf_list[i].name);
      cart_error("ART is terminating.");
    }
}


void control_parameter_list_imf(FILE *stream, const void *ptr)
{
  control_parameter_list_string(stream,imf_internal.name);
}


void config_init_imf()
{
  ControlParameterOps control_parameter_imf = { control_parameter_set_imf, control_parameter_list_imf };

  /*
  //  Default values
  */
  imf_internal.name = imf_list[1].name;
  imf_internal.f = imf_list[1].f;

  control_parameter_add2(control_parameter_imf,&imf_internal.f,"imf","imf.type","the name of the IMF. Valid names are 'Salpeter', 'Miller-Scalo', 'Chabrier', and 'Kroupa'.");

  control_parameter_add3(control_parameter_double,&imf_internal.min_mass,"imf:min-mass","imf.min_mass","am_stl","the minimum stellar mass in the IMF model.");

  control_parameter_add3(control_parameter_double,&imf_internal.max_mass,"imf:max-mass","imf.max_mass","am_stu","the maximum stellar mass in the IMF model.");
}


void config_verify_imf()
{
  /*
  //  IMF
  */
  VERIFY(imf, 1 );

  VERIFY(imf:min-mass, imf_internal.min_mass > 0.0 );

  VERIFY(imf:max-mass, imf_internal.max_mass > 0.0 );
}


double imf_f_Salpeter( double m )
{
  return pow( m, -2.35 );
}


double imf_f_MillerScalo( double m )
{
  /* Miller-Scalo (1979, ApJS 41, 513, eq. 30, Table 7) IMF */
  const double C_1 =  1.09;
  const double C_2 = -1.02;
  return exp( -C_1 * ( (log10(m) - C_2)*(log10(m) - C_2) ) ) / m;
}


double imf_f_Chabrier( double m )
{
  /* Chabrier, G. (2001, ApJ 554, 1274) */
  /* above 1Msun slope is unconstrained -- assume Salpeter*/
  const double m0_Ch = 716.4;
  const double beta_Ch = 0.25;
  const double alpha_Ch = -3.3;
  const double slope_highmass = -2.35;
  if(m < 1.0)
      return exp( -pow( m0_Ch/m, beta_Ch) ) * pow(m,alpha_Ch);
  else
      return exp( -pow( m0_Ch/1.0, beta_Ch) ) * pow(m/1.0,slope_highmass);
}


double imf_f_Kroupa( double m )
{
  /* Kroupa, P. (2007astro.ph..3124K) */
  if(m < 0.5)
    return pow(m/0.08,-1.3);
  else
    return pow(0.5/0.08,-1.3)*pow(m/0.5,-2.3);
}
  

double imf_fm( double mstar )
{
  return mstar * imf_internal.f(mstar);
}


double imf_fmz( double mstar )
{
  return mstar * imf_internal.f(mstar) * min( 0.2, max( 0.01*mstar - 0.06, 1e-20 ) );
}


#endif /* STAR_FORMATION */
