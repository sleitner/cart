#ifndef __EXT_SYNTHETIC_SPECTRA_H__
#define __EXT_SYNTHETIC_SPECTRA_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


typedef struct
{
  int Length;
  double *Vel;
  double *Flux;
}
ssSpectrum;


/*
//  Voigt profile and its integral between w1 and w2
*/
double voigt(double w, double wDopler, double wLorentz);
double voigtInt(double w1, double w2, double wDopler, double wLorentz);

#if defined(HYDRO) && defined(RADIATIVE_TRANSFER) && defined(COSMOLOGY)

ssSpectrum ssMakeSpectrum_HI1216(double pos[3], double theta, double phi, double len, double vel_pixel, int floor_level);

#endif /* HYDRO && RADIATIVE_TRANSFER && COSMOLOGY */

#endif
