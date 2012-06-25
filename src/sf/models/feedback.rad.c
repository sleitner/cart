#include "config.h"
#ifdef STAR_FORMATION

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
  double ion_rate;
}
rad_code;


void rad_setup(int level)
{
  rad_code.ion_rate = units->time/(3.0e6*constants->yr);
}


/*
//  Oldstyle hart source function
*/
float rad_luminosity_hart(int ipart)
{
  float x1, x2, dx, q;

  if(!particle_is_star(ipart)) return 0.0;

  /*
  //  The convention is different from HART
  */
  x1 = rad_code.ion_rate*(particle_t[ipart]-star_tbirth[ipart]-particle_dt[ipart]);
  if(x1 < 0.0) x1 = 0.0;
  if(x1 < 1.0e4)
    {
      /*
      //  This is a rough fit to Starburst99 evolving spectra
      */
      dx = rad_code.ion_rate*particle_dt[ipart];
      x2 = x1 + dx;
      q = rad_code.ion_rate*(0.8+x2*x2+x2*x1+x1*x1)/(1+x1*(0.8+x1*x1))/(1+x2*(0.8+x2*x2));

      return 1.4e-4*q;     
    }
  else
    {
      return 0.0;
    }
}


/*
//  New source function for Kroupa IMF
*/
float rad_luminosity_popM(int ipart)
{
  float x1, x2, dx, q, Z;

  if(!particle_is_star(ipart)) return 0.0;

  /*
  //  The convention is different from HART
  */
  x1 = rad_code.ion_rate*(particle_t[ipart]-star_tbirth[ipart]-particle_dt[ipart]);
  if(x1 < 0.0) x1 = 0.0;
  if(x1 < 1.0e4)
    {
      /*
      //  This is a rough fit to Starburst99 evolving spectra
      */
      dx = rad_code.ion_rate*particle_dt[ipart];
      x2 = x1 + dx;
      q = rad_code.ion_rate*(0.8+x2*x2+x2*x1+x1*x1)/(1+x1*(0.8+x1*x1))/(1+x2*(0.8+x2*x2));

#ifdef ENRICHMENT
      Z = (star_metallicity_II[ipart]
#ifdef ENRICHMENT_SNIa
	   + star_metallicity_Ia[ipart]
#endif /* ENRICHMENT_SNIa */
	   )/constants->Zsun;
      Z = MAX(1.0e-10,Z);
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

#endif /* STAR_FORMATION */
