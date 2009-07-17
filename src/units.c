#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "units.h"
#include "constants.h"
#include "tree.h"
#include "auxiliary.h"
#include "starformation.h"

#include "rt_solver.h"

double r0;
double t0;
double v0;
double rho0;
double den0;
double P0;
double T0;
double S0;
double H0;
double E0;
double aM0;

double AL_SD;

double Lbox = 1.0;

void init_units() {

	cosmology_init();
	cosmology_insure_consistency(1.0,0.0);

	/* H0 in s^-1 */
	H0 = 1e7 * cosmology->h / mpc;

	/* r0 in h^-1 Mpc (L_box = 1.0) */
	r0 = Lbox / (double)(num_grid);

	/* t0 in years */
	t0 = 2.0 / ( H0 * sqrt(cosmology->OmegaM) ) / year;

	/* v0 in km/s */
	v0 = 50.0 * r0 * sqrt(cosmology->OmegaM);

	/* rho0 in Msun/Mpc^3 */
	rho0 = ( 3 * H0 * H0 * cosmology->OmegaM / ( 8.0 * M_PI * G ) ) * (mpc*mpc*mpc/Msun);

	/* den0 in 1/cm^3 */
	den0 = 1.123e-5 * cosmology->Obh2;

	/* P0 in g/cm/s^2 (rho0*v0^2) */
	P0 = 4.697e-16 * cosmology->OmegaM * cosmology->OmegaM * r0 * r0 * cosmology->h * cosmology->h;

	/* T0 in K */
	T0 = 3.03e5 * r0 * r0 * wmu * cosmology->OmegaM;

	/* S0 in keV cm^2 (assumes gamma = 5/3) */
	S0 = 52.077 * pow(wmu,5./3.) * pow(cosmology->h,-4./3.) * pow(cosmology->OmegaM,1./3.) * r0 * r0;

	AL_SD = 1.6625 * (1.0 - Y_p) * (1.0 - Y_p) * hubble / sqrt(cosmology->OmegaM) / ( r0 * r0 );

	/* E0 in ergs */
	E0 = 1.38e58 * pow(r0,5.0) * cosmology->OmegaM * cosmology->OmegaM / cosmology->h;

	/* mass conversion */
	aM0 = rho0 * pow(Lbox/cosmology->h,3.0) / (double)num_root_cells;

#ifdef STARFORM
	/* call init_star_formation here so anytime units change star units change too */
	init_star_formation();
#endif /* STARFORM */

#ifdef RADIATIVE_TRANSFER
	rtSetTemUnits();
#endif

	if ( local_proc_id == MASTER_NODE ) {
		cart_debug("H0 = %e", H0 );
		cart_debug("r0 = %e", r0 );
		cart_debug("t0 = %e", t0 );
		cart_debug("v0 = %e", v0 );
		cart_debug("rho0 = %e", rho0 );
		cart_debug("den0 = %e", den0 );
		cart_debug("P0 = %e", P0 );
		cart_debug("T0 = %e", T0 );
		cart_debug("E0 = %e", E0 );
		cart_debug("aM0 = %e", aM0 );
	}
}
