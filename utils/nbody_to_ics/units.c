#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "units.h"
#include "constants.h"
#include "auxiliary.h"

int num_grid;

double aexpn;
double tl;

float cell_size[20];
float cell_size_inverse[20];
float cell_volume[20];
float cell_volume_inverse[20];

double r0;
double t0;
double v0;
double rho0;
double P0;
double T0;
double S0;
double H0;
double E0;
double aM0;

double AL_SD;

double Omega0 = 1.0;
double Omegab0 = 0.0;
double OmegaL0 = 0.0;
double hubble = 0.5;
double Lbox = 1.0;

void init_units() {
	int i;
	double a;

	/* H0 in s^-1 */
	H0 = 1e7 * hubble / mpc;

	/* r0 in h^-1 Mpc (L_box = 1.0) */
	r0 = Lbox / (double)(num_grid);

	/* t0 in years */
	t0 = 2.0 / ( H0 * sqrt( Omega0 ) ) / year;

	/* v0 in km/s */
	v0 = 50.0 * r0 * sqrt(Omega0);

	/* rho0 in Msun/Mpc^3 */
	rho0 = (3.0*H0*H0*Omega0 / ( 8.0 * M_PI * G ))*(mpc*mpc*mpc/Msun);

	/* P0 in g/cm/s^2 (rho0*v0^2) */
	P0 = 4.697e-16 * Omega0 * Omega0 * r0 * r0 * hubble * hubble;

	/* T0 in K */
	T0 = 3.03e5 * r0*r0 * wmu * Omega0;

	/* S0 in keV cm^2 (assumes gamma = 5/3) */
	S0 = 52.077 * pow( wmu, 5./3.) * pow( hubble, -4./3. )*pow(Omega0,1./3.)*r0*r0;

	AL_SD = 1.6625 * (1.0 - Y_p)*(1.0 - Y_p) * hubble / sqrt(Omega0) / (r0*r0);

	/* E0 in ergs */
	E0 = 1.38e58 * pow(r0,5.0) * Omega0*Omega0 / hubble;

	/* mass conversion */
	aM0 = rho0 * pow( Lbox / hubble, 3.0 ) / ((double)num_grid*(double)num_grid*(double)num_grid);

        for ( i = 0; i < 20; i++ ) {
                cell_size_inverse[i] = (float)(1<<i);
                cell_size[i] = 1.0 / cell_size_inverse[i];
                cell_volume_inverse[i] = (float)(1 << (nDim*i));
                cell_volume[i] = 1.0 / cell_volume_inverse[i];
	}

	tl = a2b( aexpn );
}

double f_a2b( double a ) {
	return 0.5*sqrt(Omega0) / ( a*a*a*sqrt(Omega0/(a*a*a) + OmegaL0 + (1.0-Omega0-OmegaL0)/(a*a) ) );
}

/* convert expansion factor a into hydro timestep var b */
double a2b( double a ) {

	cart_assert( a >= 0 );
	cart_assert( Omega0 > 0.0 || OmegaL0 > 0.0 );

	if ( Omega0 == 1.0 && OmegaL0 == 0.0 ) {
		return (1.0 - 1.0/sqrt(a));
	} else {
		return integrate( f_a2b, 1.0, a, 1e-9, 1e-9 );
	}
}

double b0;
double f_b2a( double a ) {
	return a2b(a) - b0;
}

double b2a( double b ) {
	double a;

	cart_assert ( Omega0 > 0.0 || OmegaL0 > 0.0 );

	if ( Omega0 == 1.0 && OmegaL0 == 0.0 ) {
		a = (1.0/(1.0 - b))*(1.0/(1.0 - b));
	} else {
		b0 = b;
		a = root_finder( f_b2a, 1e-4, 1.1, 1e-9, 1e-9 );
	}

	return a; 
}

double f_age( double a ) {
	double fact;

	fact = Omega0 + OmegaL0*a*a*a + (1.0 - Omega0 - OmegaL0)*a;

	if ( fact > 0.0 ) {
		return sqrt(a/fact);
	} else {
		return 0.0;
	}
}

double age( double a ) {
	return integrate( f_age, 0.0, a, 1e-9, 1e-9 ) / ( 1e7 * hubble * gyr/mpc );
}

double age_t( double t ) {
	double a;
	double f;

	a = b2a(t);
	f = Omega0 / ( 1.0 - Omega0 ) / (a*a*a);
	return 9.779/hubble * 2.0/3.0/
		sqrt(1.0-Omega0)*log((1.0+sqrt(1.0+f))/sqrt(f));
}
