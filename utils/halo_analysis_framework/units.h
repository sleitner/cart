#ifndef __UNITS_H__
#define __UNITS_H__

#include "constants.h"

extern int num_grid;
extern double aexpn;
extern double tl;

extern double r0;
extern double t0;
extern double v0;
extern double rho0;
extern double rho0cgs;
extern double P0;
extern double T0;
extern double S0;
extern double H0;
extern double E0;
extern double aM0;

extern double AL_SD;

extern double Omega0;
extern double Omegab0;
extern double OmegaL0;
extern double hubble;
extern double Lbox;

extern double cell_size[20];
extern double cell_size_inverse[20];
extern double cell_volume[20];
extern double cell_volume_inverse[20];

void init_units();
double a2b( double a );
double b2a( double b );
double age( double a );
double age_t( double t );

#endif
