#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "auxiliary.h"
#include "analysis.h"
#include "units.h"
#include "constants.h"
#include "cooling.h"
#include "io.h"
#include "sfc.h"

double compute_distance_periodic( double pos1[nDim], double pos2[nDim] ) {
	int d;
	double dx, r = 0.0;

	for ( d = 0; d < nDim; d++ ) {
		dx = fabs( pos1[d] - pos2[d] );

		if ( dx > (double)(num_grid/2) ) {
			dx -= num_grid;
		}

		r += dx*dx;
	}

	return sqrt(r);
}

double fact_nH;
double Tfact, vfact, rfact, Sfact;
double szfact, Pfact;
double nfact;

void cell_callback( halo_struct *halo, cell_struct *cell ) {
	int i, j;
	int level;
	double Tcell, Pcell, Scell, Ycell;
	double rhogi;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double Tecell, Tecell_kev, Pecell, Yecell;
	double ne, loglambda, tei;
#endif

	level = cell->level;

	rhogi = 1.0 / cell->gas_density;

	/* cell properties (K, ergs cm^-2, keV cm^-2, Mpc^2) */
	Tcell = Tfact * cell->gas_internal_energy * rhogi;
	Pcell = Pfact * cell->gas_pressure;
	Scell = Sfact * cell->gas_internal_energy*pow(rhogi,gamma);
	Ycell = szfact * cell->gas_internal_energy;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	Pecell = Pfact * (2./3.)*cell->electron_internal_energy;
	Tecell = Tfact * (wmu_e/wmu) * cell->electron_internal_energy * rhogi;
	Yecell = szfact * (wmu_e/wmu) * cell->electron_internal_energy;
	ne = nfact*cell->gas_density;
	tei = 6.3e8*pow(Tecell/1e7,1.5)/(ne/1e-5)/(loglambda/40.0);
#endif

}

void particle_callback( halo_struct *halo, particle_struct *particle ) {
	profile_struct *profile;
}

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos, halo_list *subhalos ) {
	int i, j, k;
	char filename1[256];
	char filename2[256];
	char filename3[256];

	/* set up conversion constants */
#ifdef HYDRO
	Tfact = T0 * ( gamma - 1.0 ) / ( aexpn*aexpn );
	Pfact = P0/(aexpn*aexpn*aexpn*aexpn*aexpn);
	Sfact = S0; /* no trend with aexp if gamma = 5/3 */
	szfact = 3.4383e-15 * ( gamma - 1.0 ) * T0 * r0 * Omega0 * hubble / (aexpn*aexpn*aexpn*aexpn);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	nfact = 1.12e-5*hubble*hubble*Omega0 / wmu_e / (aexpn*aexpn*aexpn);
#endif 
#endif

	rfact = 1000.0 * r0 * aexpn / hubble; /* code units -> proper kpc */
	vfact = v0 / aexpn;

#ifdef PARTICLES
	/* open particle file */
	sprintf( filename1, "%s/PMcrda%06.4f.DAT", output_directory, aexpn );
	sprintf( filename2, "%s/PMcrs_indexed_a%06.4f.DAT", output_directory, aexpn );
	sprintf( filename3, "%s/stars_indexed_a%06.4f.dat", output_directory, aexpn );

	read_indexed_particles( filename1, filename2, filename3, halos, subhalos, particle_callback );

	cart_debug("done assigning particles");
#endif /* PARTICLES */

#ifdef HYDRO
	sprintf(filename1, "%s/%s_a%6.4f.grid", output_directory, jobname, aexpn );
	read_indexed_grid( filename1, halos, subhalos, cell_callback );

	cart_debug("done assigning gas to grid"); fflush(stdout);
#endif /* HYDRO */

}
