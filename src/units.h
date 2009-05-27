#ifndef __UNITS_H__
#define __UNITS_H__

#include "constants.h"
#include "cosmology.h"

extern double r0;
extern double t0;
extern double v0;
extern double rho0;
extern double den0;
extern double P0;
extern double T0;
extern double S0;
extern double H0;
extern double E0;
extern double aM0;

extern double AL_SD;
extern double Lbox;

void init_units();

/* NG: legacy names for backward compatibility */
#define Omega0      cosmology->OmegaM
#define Omegab0     cosmology->OmegaB
#define OmegaL0     cosmology->OmegaL
#define hubble      cosmology->h

#endif
