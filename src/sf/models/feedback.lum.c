#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "particle.h"
#include "starformation.h"
#include "units.h"


/*
//  Radiative feedback
*/

struct
{
  double ion_time_scale;
}
lum = { 3.0e6 };


void lum_config_init()
{
  control_parameter_add2(control_parameter_time,&lum.ion_time_scale,"lum:ion-time-scale","lum.ion_time_scale","time-scale for the evolution of the ionizing luminosity from young stars.");
}


void lum_config_verify()
{
  /*
  //  ionizing luminosity
  */
  VERIFY(lum:ion-time-scale, lum.ion_time_scale > 0.0 );
}


struct
{
  double ion_rate;
}
lum_code;


void lum_setup(int level)
{
  lum_code.ion_rate = units->time/(lum.ion_time_scale*constants->yr);
}


float lum_ionizing_luminosity(int ipart)
{
  float x1, x2, dx, q;

  if(!particle_is_star(ipart)) return 0.0;

  /*
  //  The convention is different from HART
  */
  x1 = lum_code.ion_rate*(particle_t[ipart]-star_tbirth[ipart]);
  if(x1 < 0.0) x1 = 0.0;
  if(x1 < 1.0e4)
    {
      /*
      //  This is a rough fit to Starburst99 evolving spectra
      */
      dx = lum_code.ion_rate*particle_dt[ipart];
      if(dx > 1.0e-5)
        {
          x2 = x1 + dx;
          x1 *= (0.8+x1*x1);
          x2 *= (0.8+x2*x2);
          q = (x2-x1)/(1+x1)/(1+x2)/particle_dt[ipart];
        }
      else
        {
          x2 = x1*(0.8+x1*x1);
          q = (0.8+3*x1*x1)/(1+x2)/(1+x2)*lum_code.ion_rate;
        }
      return 1.4e-4*q;     
    }
  else
    {
      return 0.0;
    }
}


float lum2012_ionizing_luminosity(int ipart)
{
  float x1, x2, dx, q, Z;

  if(!particle_is_star(ipart)) return 0.0;

  /*
  //  The convention is different from HART
  */
  x1 = lum_code.ion_rate*(particle_t[ipart]-star_tbirth[ipart]);
  if(x1 < 0.0) x1 = 0.0;
  if(x1 < 1.0e4)
    {
      /*
      //  This is a rough fit to Starburst99 evolving spectra
      */
      dx = lum_code.ion_rate*particle_dt[ipart];
      if(dx > 1.0e-5)
        {
          x2 = x1 + dx;
          x1 *= (0.8+x1*x1);
          x2 *= (0.8+x2*x2);
          q = (x2-x1)/(1+x1)/(1+x2)/particle_dt[ipart];
        }
      else
        {
          x2 = x1*(0.8+x1*x1);
          q = (0.8+3*x1*x1)/(1+x2)/(1+x2)*lum_code.ion_rate;
        }
#ifdef ENRICHMENT
      Z = (star_metallicity_II[ipart]
#ifdef ENRICHMENT_SNIa
	   + star_metallicity_Ia[ipart]
#endif /* ENRICHMENT_SNIa */
	   )/constants->Zsun;
      Z = max(1.0e-10,Z);
#else  /* ENRICHMENT */
      Z = 0.1;
#endif /* ENRICHMENT */
      return 1.04e-4/powf(Z,0.1f)/(1.0+0.27*Z)*q;     
    }
  else
    {
      return 0.0;
    }
}

#endif /* STARFORM */
