#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "auxiliary.h"
#include "analysis.h"
#include "analysis_xray.h"
#include "units.h"
#include "constants.h"
#include "cooling.h"
#include "io.h"
#include "sfc.h"

int num_radii = 4;
char *radii_label[] = { 
	"rout", /* rout MUST be first radius */
	"r500c",
	"r200c",
	"r180m" };

float radii_delta[] = {	
	0.0,
	500.0,
	200.0,
	180.0 };

char radii_units[] = {	
	'r',
	'c',
	'c',
	'm' };

double rr[max_bins], rl[max_bins], bin_volume[max_bins], bin_volume_cumulative[max_bins];

double log_interpolate( double *binned_var, int bin, double rlout, double rri, double rll ) {
	double aM1, aM2, ah, bh;

	if ( binned_var[bin-1] > 0.0 ) {
		aM1 = log10( binned_var[bin-1] );
	} else {
		aM1 = -15.0;
	}
	if ( binned_var[bin] > 0.0 ) {
		aM2 = log10( binned_var[bin] );
	} else {
		aM2 = -15.0;
	}
	ah = ( aM2 - aM1 ) * rri;
	bh = aM1 - ah * rll;
	return pow( 10.0, ah*rlout + bh );
}

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
double dtfact, dEfact, szfact, Pfact;
double tnewstar;
double nfact;

int compute_radial_bin( double r, double rlmin, double drl ) {
	int bin;
	double lgr = log10(r);

	if ( lgr < rlmin ) {
		bin = 0;
	} else {
		bin = (int)((lgr - rlmin)/drl) + 1;
	}

	return bin;
}

double corner_offset[8][nDim] = {
	{ -0.5, -0.5, -0.5 }, {  0.5, -0.5, -0.5 },
	{ -0.5,  0.5, -0.5 }, {  0.5,  0.5, -0.5 },
	{ -0.5, -0.5,  0.5 }, {  0.5, -0.5,  0.5 },
	{ -0.5,  0.5,  0.5 }, {  0.5,  0.5,  0.5 }
};

void cell_callback( halo_struct *halo, cell_struct *cell ) {
	int i, j;
	int bin, bin2;
	int differ;
	int level;
	double cell_mass;
	double Tcell, Tcell_kev;
	double Pcell, Scell;
	double tcool, dEcell;
	double Ycell;
	double rhogi, rhogl;
	double Zdum, Zldum;
	double cell_vx, cell_vy, cell_vz;

	double Fline, Fcont, Tline, Tcont, xx, f_line;
	double Tcont1, Tcont2, avgE;
	double xray_Tcont1_cell, xray_Tcont2_cell;
	double xray_cT, xray_lambda, xray_fT, xray_w;
	double xray_Tcont1, xray_Tcont2;
	double xray_Fcont_cell, xray_Fline_cell, xray_avgE_cell;

	double r;
	double pos[nDim];
	double mass_fraction, volume_fraction;
	int bin_volume_fraction[max_bins];
	int num_points;
	int min_bin;
	double min_bin_volume;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double Tecell, Tecell_kev, Pecell, Yecell;
	double ne, loglambda, tei;

	double eFline, eFcont, eTline, eTcont, exx, ef_line;
	double eTcont1, eTcont2, eavgE;
	double exray_Tcont1_cell, exray_Tcont2_cell;
	double exray_cT, exray_lambda, exray_fT, exray_w;
	double exray_Tcont1, exray_Tcont2;
	double exray_Fcont_cell, exray_Fline_cell, exray_avgE_cell;
#endif
	profile_struct *profile;

	profile = (profile_struct *)halo->data;

	level = cell->level;

	r = compute_distance_periodic( halo->pos, cell->pos );
	bin = compute_radial_bin( r, profile->rlmin, profile->drl );
	min_bin = bin;

	differ = 0;
	for ( i = 0; i < 8; i++ ) {
		for ( j = 0; j < nDim; j++ ) {
			pos[j] = cell->pos[j] + corner_offset[i][j]*cell_size[level];
		}

		r = compute_distance_periodic( halo->pos, pos );
		bin2 = compute_radial_bin( r, profile->rlmin, profile->drl );

		if ( bin2 < profile->num_bins && bin != bin2 ) {
			differ = 1;
		}

		min_bin = min( min_bin, bin2 );
	}

	if ( !differ && min_bin >= profile->num_bins ) {
		return;
	}

	rhogi = 1.0 / cell->gas_density;

	/* grab cell properties */
	Tcell = Tfact * cell->gas_internal_energy * rhogi;
	Tcell_kev = Tcell * keV;

	Pcell = Pfact * cell->gas_pressure;
	Scell = Sfact * cell->gas_internal_energy*pow(rhogi,gamma);
	Ycell = szfact * cell->gas_internal_energy;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	Pecell = Pfact * (2./3.)*cell->electron_internal_energy;
	Tecell = Tfact * (wmu_e/wmu) * cell->electron_internal_energy * rhogi;
	Tecell_kev = Tecell * keV;
	Yecell = szfact * (wmu_e/wmu) * cell->electron_internal_energy;
	ne = nfact*cell->gas_density;
	loglambda = max( 30.0, 37.8 + log(Tecell/1e7) - 0.5*log(ne/1e-5) );
	tei = 6.3e8*pow(Tecell/1e7,1.5)/(ne/1e-5)/(loglambda/40.0);
#endif

#ifdef COOLING
#ifdef CLOUDY_COOLING
	/* take code density -> log10(n_H [cm^-3]) */
	rhogl = log10(cell->gas_density) + fact_nH;
#endif /* CLOUDY_COOLING */

#ifdef METALCOOLING
	Zdum = cell->metallicity_II+cell->metallicity_Ia;
	Zdum = max( Zdum, 1e-10 );
	Zdum *= rhogi / Zsolar;
	Zldum = log10(Zdum);
#else
	Zdum = 0.0;
	Zldum = 0.0;
#endif /* METALCOOLING*/

	dEcell = cooling_rate( rhogl, Tcell*1e-4, Zldum ) * 
		cell->gas_density*cell->gas_density * aexpn;
	tcool = cell->gas_internal_energy / dEcell / dtfact;
	dEcell /= dEfact;
#else
	Zdum = 0.3;
	Zldum = log10(Zdum);
#endif /* COOLING */

	cell_mass = cell->gas_density*cell_volume[level];

	cell_vx = cell->momentum[0]*rhogi;
	cell_vy = cell->momentum[1]*rhogi;
	cell_vz = cell->momentum[2]*rhogi;

#ifdef ANALYSIS_XRAY
	/* Alexey's weighting (Vikhlinin 2006) */
	xray_calibration( Tcell_kev, &xray_cT, &xray_lambda, &xray_fT );

	/* eq 6 */
	xray_w = xray_cT * cell->gas_density*cell->gas_density*pow(Tcell_kev,-alpha_xray);

	/* numerator & denominator of eq 4, computes <T>_cont */
	xray_Tcont1 = xray_w*Tcell_kev*cell_volume[level];
	xray_Tcont2 = xray_w*cell_volume[level];

	/* eq 9 */
	xray_Fcont_cell = xray_cT * cell->gas_density *
		cell->gas_density*cell_volume[level];
	cart_assert( xray_cT >= 0.0 );
	cart_assert( xray_Fcont_cell >= 0.0 );

	/* eq 11 */
	xray_Fline_cell = xray_lambda * Zdum *
		cell->gas_density*cell->gas_density *
		cell_volume[level];
	xray_avgE_cell = xray_fT * xray_Fline_cell;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	xray_calibration( Tecell_kev, &exray_cT, &exray_lambda, &exray_fT );
	exray_w = exray_cT * cell->gas_density*cell->gas_density*pow(Tecell_kev,-alpha_xray);
	exray_Tcont1 = exray_w*Tecell_kev*cell_volume[level];
	exray_Tcont2 = exray_w*cell_volume[level];
	exray_Fcont_cell = exray_cT * cell->gas_density *
		cell->gas_density*cell_volume[level];
	exray_Fline_cell = exray_lambda * Zdum *
		cell->gas_density*cell->gas_density *
		cell_volume[level];
	exray_avgE_cell = exray_fT * exray_Fline_cell;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */

	if ( !differ ) {
		profile->bin_gas_mass[bin] += cell_mass;

		if ( !cell->subhalo_flag ) {
			profile->bin_excluded_gas_mass[bin] += cell_mass;
			profile->bin_excluded_volume[bin] += cell_volume[level];

			if ( Tcell <= Tcold ) {
				profile->bin_cold_gas_mass[bin] += cell_mass;
			}

			profile->bin_sz_flux[bin] += Ycell*cell_volume[level];
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			profile->bin_electron_sz_flux[bin] += Yecell*cell_volume[level];
#endif

#ifdef MASS_WEIGHTED
			profile->bin_massweighted_pressure[bin] += Pcell*cell_mass;
			profile->bin_massweighted_entropy[bin] += Scell*cell_mass;
			profile->bin_massweighted_temperature[bin] += Tcell*cell_mass;
#endif /* MASS_WEIGHTED */

			profile->bin_gas_pressure[bin] += Pcell*cell_volume[level];
			profile->bin_gas_entropy[bin] += Scell*cell_volume[level];
			profile->bin_gas_temperature[bin] += Tcell*cell_volume[level];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			profile->bin_electron_pressure[bin] += Pecell*cell_volume[level];
			profile->bin_electron_temperature[bin] += Tecell*cell_volume[level];
			profile->bin_electron_density[bin] += ne*cell_volume[level];
			profile->bin_electron_tei[bin] += tei*cell_volume[level];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef COOLING
			profile->bin_gas_coolingrate[bin] += dEcell*cell_volume[level];
			profile->bin_gas_tcool[bin] += tcool;
#endif /* COOLING */

#ifdef ENRICH
			profile->bin_gas_metallicity_II[bin] += cell->metallicity_II*rhogi*cell_mass;
			profile->bin_gas_metallicity_Ia[bin] += cell->metallicity_Ia*rhogi*cell_mass;
#endif /* ENRICH */

			profile->bin_gas_velocity[0][bin] += cell_vx*cell_mass;
			profile->bin_gas_velocity[1][bin] += cell_vy*cell_mass;
			profile->bin_gas_velocity[2][bin] += cell_vz*cell_mass;

			profile->bin_gas_vrms[bin] += (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*cell_mass;

#ifdef ANALYSIS_XRAY
			/* X-ray quantities */
			if ( Tcell_kev > 0.086 ) {
				profile->bin_xray_Fcont[bin] += xray_Fcont_cell;
				profile->bin_xray_Fline[bin] += xray_Fline_cell;
				profile->bin_xray_avgE[bin] += xray_avgE_cell;
				profile->bin_xray_Tcont1[bin] += xray_Tcont1;
				profile->bin_xray_Tcont2[bin] += xray_Tcont2;
			}
#ifdef ELECTRON_ION_NONEQUILIBRIUM 
			if ( Tecell_kev > 0.086 ) {
				profile->bin_exray_Fcont[bin] += exray_Fcont_cell;
				profile->bin_exray_Fline[bin] += exray_Fline_cell;
				profile->bin_exray_avgE[bin] += exray_avgE_cell;
				profile->bin_exray_Tcont1[bin] += exray_Tcont1;
				profile->bin_exray_Tcont2[bin] += exray_Tcont2;
			}
#endif /* ELECGTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */
		}
	} else {
		/* special case of halo center lying inside cell */
		if ( 	!cell->subhalo_flag && 
				halo->pos[0] > cell->pos[0] - 0.5*cell_size[level] && 
				halo->pos[0] < cell->pos[0] + 0.5*cell_size[level] &&
				halo->pos[1] > cell->pos[1] - 0.5*cell_size[level] &&
				halo->pos[1] < cell->pos[1] + 0.5*cell_size[level] &&
				halo->pos[2] > cell->pos[2] - 0.5*cell_size[level] &&
				halo->pos[2] < cell->pos[2] + 0.5*cell_size[level]  ) {
			/* find largest bin entirely contained within cell */
			r = halo->pos[0] - cell->pos[0] + 0.5*cell_size[level];	
			r = min( r, cell->pos[0] + 0.5*cell_size[level] - halo->pos[0] );
			r = min( r, halo->pos[1] - cell->pos[1] + 0.5*cell_size[level] );
			r = min( r, cell->pos[1] + 0.5*cell_size[level] - halo->pos[1] );
			r = min( r, halo->pos[2] - cell->pos[2] + 0.5*cell_size[level] );
			r = min( r, cell->pos[2] + 0.5*cell_size[level] - halo->pos[2] );

			min_bin = compute_radial_bin( r, profile->rlmin, profile->drl );

			for ( bin = 0; bin < min_bin; bin++ ) {
				volume_fraction = bin_volume[bin];
				mass_fraction = bin_volume[bin]*cell->gas_density;

				profile->bin_gas_mass[bin] = mass_fraction;
				profile->bin_excluded_gas_mass[bin] = mass_fraction;
				profile->bin_volume[bin] = volume_fraction;
				profile->bin_excluded_volume[bin] = volume_fraction;

				if ( Tcell <= Tcold ) {
					profile->bin_cold_gas_mass[bin] = mass_fraction;
				}

				profile->bin_sz_flux[bin] = Ycell*volume_fraction;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				profile->bin_electron_sz_flux[bin] = Yecell*volume_fraction;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef MASS_WEIGHTED
				profile->bin_massweighted_pressure[bin] = Pcell*mass_fraction;
				profile->bin_massweighted_entropy[bin] = Scell*mass_fraction;
				profile->bin_massweighted_temperature[bin] = Tcell*mass_fraction;
#endif /* MASS_WEIGHTED */

				profile->bin_gas_pressure[bin] = Pcell*volume_fraction;
				profile->bin_gas_entropy[bin] = Scell*volume_fraction;
				profile->bin_gas_temperature[bin] = Tcell*volume_fraction;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
				profile->bin_electron_pressure[bin] = Pecell*volume_fraction;
				profile->bin_electron_temperature[bin] = Tecell*volume_fraction;
				profile->bin_electron_density[bin] = ne*volume_fraction;
				profile->bin_electron_tei[bin] = tei*volume_fraction;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef COOLING
				profile->bin_gas_coolingrate[bin] = dEcell*volume_fraction;
				profile->bin_gas_tcool[bin] = tcool*mass_fraction;
#endif /* COOLING */

#ifdef ENRICH
				profile->bin_gas_metallicity_II[bin] = cell->metallicity_II*rhogi*mass_fraction;
				profile->bin_gas_metallicity_Ia[bin] = cell->metallicity_Ia*rhogi*mass_fraction;
#endif /* ENRICH */

				profile->bin_gas_velocity[0][bin] = cell_vx*mass_fraction;
				profile->bin_gas_velocity[1][bin] = cell_vy*mass_fraction;
				profile->bin_gas_velocity[2][bin] = cell_vz*mass_fraction;

				profile->bin_gas_vrms[bin] = (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*mass_fraction;

#ifdef ANALYSIS_XRAY
				/* X-ray quantities */
				profile->bin_xray_Fcont[bin] = xray_Fcont_cell*volume_fraction;
				profile->bin_xray_Fline[bin] = xray_Fline_cell*volume_fraction;
				profile->bin_xray_avgE[bin] = xray_avgE_cell*volume_fraction;
				profile->bin_xray_Tcont1[bin] = xray_Tcont1*volume_fraction;
				profile->bin_xray_Tcont2[bin] = xray_Tcont2*volume_fraction;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				profile->bin_exray_Fcont[bin] = exray_Fcont_cell*volume_fraction;
				profile->bin_exray_Fline[bin] = exray_Fline_cell*volume_fraction;
				profile->bin_exray_avgE[bin] = exray_avgE_cell*volume_fraction;
				profile->bin_exray_Tcont1[bin] = exray_Tcont1*volume_fraction;
				profile->bin_exray_Tcont2[bin] = exray_Tcont2*volume_fraction;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */
			}
		}

		for ( bin = 0; bin < profile->num_bins; bin++ ) {
			bin_volume_fraction[bin] = 0;
		}

		/* assure reasonable coverage if cell is large fraction of minimum bin */
		min_bin_volume = 4.*M_PI/3.*pow(10.,3.0*(profile->rlmin+(double)min_bin*profile->drl));
		if ( cell_volume[level] > 0.1*min_bin_volume ) {
			num_points = 10*points_per_cell*cell_volume[level]/min_bin_volume;
			cart_assert( num_points >= points_per_cell );
		} else {
			num_points = points_per_cell;
		}

		/* throw random points */
		for ( i = 0; i < num_points; i++ ) {
			for ( j = 0; j < nDim; j++ ) {
				pos[j] = cell->pos[j] + cell_size[level]*(cart_rand()-0.5);
			}

			r = compute_distance_periodic( halo->pos, pos );
			bin = compute_radial_bin( r, profile->rlmin, profile->drl );

			if ( bin < profile->num_bins ) {	
				bin_volume_fraction[bin]++;
			}
		}

		for ( bin = min_bin; bin < profile->num_bins; bin++ ){ 
			if ( bin_volume_fraction[bin] > 0 ) {
				volume_fraction = (double)bin_volume_fraction[bin] /
					(double)num_points;
				cart_assert( volume_fraction > 0.0 && volume_fraction <= 1.0 );

				mass_fraction = cell_mass*volume_fraction;
				volume_fraction *= cell_volume[level];

				profile->bin_gas_mass[bin] += mass_fraction;
				profile->bin_volume[bin] += volume_fraction;

				if ( !cell->subhalo_flag ) {
					profile->bin_excluded_gas_mass[bin] += mass_fraction;
					profile->bin_excluded_volume[bin] += volume_fraction;

					if ( Tcell <= Tcold ) {
						profile->bin_cold_gas_mass[bin] += mass_fraction;
					}

					profile->bin_sz_flux[bin] += Ycell*volume_fraction;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
					profile->bin_electron_sz_flux[bin] += Yecell*volume_fraction;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef MASS_WEIGHTED
					profile->bin_massweighted_pressure[bin] += Pcell*mass_fraction;
					profile->bin_massweighted_entropy[bin] += Scell*mass_fraction;
					profile->bin_massweighted_temperature[bin] += Tcell*mass_fraction;
#endif /* MASS_WEIGHTED */

					profile->bin_gas_pressure[bin] += Pcell*volume_fraction;
					profile->bin_gas_entropy[bin] += Scell*volume_fraction;
					profile->bin_gas_temperature[bin] += Tcell*volume_fraction;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					profile->bin_electron_pressure[bin] += Pecell*volume_fraction;
					profile->bin_electron_temperature[bin] += Tecell*volume_fraction;
					profile->bin_electron_density[bin] += ne*volume_fraction;
					profile->bin_electron_tei[bin] += tei*volume_fraction;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef COOLING
					profile->bin_gas_coolingrate[bin] += dEcell*volume_fraction;
					profile->bin_gas_tcool[bin] += tcool*mass_fraction;
#endif /* COOLING */

#ifdef ENRICH
					profile->bin_gas_metallicity_II[bin] += cell->metallicity_II*rhogi*mass_fraction;
					profile->bin_gas_metallicity_Ia[bin] += cell->metallicity_Ia*rhogi*mass_fraction;
#endif /* ENRICH */

					profile->bin_gas_velocity[0][bin] += cell_vx*mass_fraction;
					profile->bin_gas_velocity[1][bin] += cell_vy*mass_fraction;
					profile->bin_gas_velocity[2][bin] += cell_vz*mass_fraction;

					profile->bin_gas_vrms[bin] += (cell_vx*cell_vx + cell_vy*cell_vy + cell_vz*cell_vz)*mass_fraction;

#ifdef ANALYSIS_XRAY
					/* X-ray quantities */
					profile->bin_xray_Fcont[bin] += xray_Fcont_cell*volume_fraction;
					profile->bin_xray_Fline[bin] += xray_Fline_cell*volume_fraction;
					profile->bin_xray_avgE[bin] += xray_avgE_cell*volume_fraction;
					profile->bin_xray_Tcont1[bin] += xray_Tcont1*volume_fraction;
					profile->bin_xray_Tcont2[bin] += xray_Tcont2*volume_fraction;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
					profile->bin_exray_Fcont[bin] += exray_Fcont_cell*volume_fraction;
					profile->bin_exray_Fline[bin] += exray_Fline_cell*volume_fraction;
					profile->bin_exray_avgE[bin] += exray_avgE_cell*volume_fraction;
					profile->bin_exray_Tcont1[bin] += exray_Tcont1*volume_fraction;
					profile->bin_exray_Tcont2[bin] += exray_Tcont2*volume_fraction;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */
				}
			}
		}
	}
}

void particle_callback( halo_struct *halo, particle_struct *particle ) {
	int bin;
	double r;
	profile_struct *profile;

	profile = (profile_struct *)halo->data;

	r = compute_distance_periodic( halo->pos, particle->x );
	bin = compute_radial_bin( r, profile->rlmin, profile->drl );

	if ( bin >= 0 && bin < profile->num_bins ) {

#ifdef STARFORM
		if ( particle->is_star ) {
			profile->bin_star_num[bin]++;
			profile->bin_star_mass[bin] += particle->mass;

			profile->bin_star_momentum[0][bin] += particle->v[0]*particle->mass;
			profile->bin_star_momentum[1][bin] += particle->v[1]*particle->mass;
			profile->bin_star_momentum[2][bin] += particle->v[2]*particle->mass;

			profile->bin_star_vrms[bin] += particle->mass * (
					particle->v[0]*particle->v[0] +
					particle->v[1]*particle->v[1] +
					particle->v[2]*particle->v[2] );

			profile->bin_star_age[bin] += particle->star_tbirth*particle->mass;
			profile->bin_star_metallicity_II[bin] += particle->star_metallicity_II*particle->mass;
			profile->bin_star_metallicity_Ia[bin] += particle->star_metallicity_Ia*particle->mass;

			if ( tl - particle->star_tbirth <= tnewstar ) {
				profile->bin_new_star_num[bin]++;
				profile->bin_new_star_mass[bin] += particle->mass;

				profile->bin_new_star_momentum[0][bin] += particle->v[0]*particle->mass;
				profile->bin_new_star_momentum[1][bin] += particle->v[1]*particle->mass;
				profile->bin_new_star_momentum[2][bin] += particle->v[2]*particle->mass;

				profile->bin_new_star_vrms[bin] += particle->mass * (
						particle->v[0]*particle->v[0] +
						particle->v[1]*particle->v[1] +
						particle->v[2]*particle->v[2] );
				profile->bin_new_star_metallicity_II[bin] += particle->star_metallicity_II*particle->mass;
				profile->bin_new_star_metallicity_Ia[bin] += particle->star_metallicity_Ia*particle->mass;
			}
		} else {
			profile->bin_dark_num[bin]++;
			profile->bin_dark_mass[bin] += particle->mass;

			profile->bin_dark_momentum[0][bin] += particle->v[0]*particle->mass;
			profile->bin_dark_momentum[1][bin] += particle->v[1]*particle->mass;
			profile->bin_dark_momentum[2][bin] += particle->v[2]*particle->mass;

			profile->bin_dark_vrms[bin] += particle->mass * (
					particle->v[0]*particle->v[0] +
					particle->v[1]*particle->v[1] +
					particle->v[2]*particle->v[2] );

		}
#else
		profile->bin_dark_num[bin]++;
		profile->bin_dark_mass[bin] += particle->mass;
		profile->bin_dark_momentum[0][bin] += particle->v[0]*particle->mass;
		profile->bin_dark_momentum[1][bin] += particle->v[1]*particle->mass;
		profile->bin_dark_momentum[2][bin] += particle->v[2]*particle->mass;

		profile->bin_dark_vrms[bin] += particle->mass * (
				particle->v[0]*particle->v[0] +
				particle->v[1]*particle->v[1] +
				particle->v[2]*particle->v[2] );
#endif /* STARFORM */
	}

}

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos, halo_list *subhalos ) {
	int i, j, k;
	int ihalo, icell, ipart;
	int ix, iy, iz;
	int index;
	int level;
	int bin;
	double r, rrl, rll, rri, rout, rlout;
	int max_bin, min_bin;
	int num_bins;
	double rlmin, rlmax, drl, dlout;
	double rlvirmin, rlvIrmax, drlvir;
	double rvir_bin_volume[max_bins];
	double rmid[max_bins];
	double inverse_mass, vmean;
	int irvir, irflag;
	double rvir, rdout, rvdout, rmass;
	double dbi1, dbi2, dlbi1, dlbi2;
	double aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_dark, aM_baryons, aM_total;
	double Ysz, Ysz_total, Tx;
	double total_dark_mass;
	double total_star_mass;
	double total_new_star_mass;
	double total_gas_mass;
	double total_cold_gas_mass;
	double avg_gas_metallicity_II;
	double avg_gas_metallicity_Ia;
	double avg_star_metallicity_II;
	double avg_star_metallicity_Ia;
	double avg_new_star_metallicity_II;
	double avg_new_star_metallicity_Ia;
	double avg_star_age;
	double avg_gas_mass;
	double avg_new_star_mass;
	double avg_star_mass;
	double xz, omega, a3, E2;
	double virial_overdensity;
	double overdensity;
	double vmax, rmax;
	double rcirc, vcirc_gas, vcirc_dm, vcirc_stars, vcirc_total;
	double age_star;

	profile_struct *profile;

	double fact_nH;
	double Fline, Fcont, Tline, Tcont, xx, f_line;
	double Tcont1, Tcont2, avgE;

	double radii_overdensity[num_radii];
	int iflag_r[num_radii];
	double delta_r[num_radii];
	double aM_gas_r[num_radii];
	double aM_cold_gas_r[num_radii];
	double aM_stars_r[num_radii];
	double aM_new_stars_r[num_radii];
	double aM_dark_r[num_radii];
	double aM_baryons_r[num_radii];
	double aM_total_r[num_radii];
	double Ysz_r[num_radii];
	double Tx_r[num_radii];
	double avg_gas_metallicity_II_r[num_radii];
	double avg_gas_metallicity_Ia_r[num_radii];
	double avg_star_metallicity_II_r[num_radii];
	double avg_star_metallicity_Ia_r[num_radii];
	double avg_new_star_metallicity_II_r[num_radii];
	double avg_new_star_metallicity_Ia_r[num_radii];
	double avg_star_age_r[num_radii];
	double avg_gas_mass_r[num_radii];
	double avg_new_star_mass_r[num_radii];
	double avg_star_mass_r[num_radii];

	double bin_total_dark_mass[max_bins];
	double bin_total_star_mass[max_bins];
	double bin_total_new_star_mass[max_bins];
	double bin_total_gas_mass[max_bins];
	double bin_total_cold_gas_mass[max_bins];
	double bin_total_density[max_bins];
	double bin_total_sz_flux[max_bins];

	double bin_total_xray_Fcont[max_bins];
	double bin_total_xray_Fline[max_bins];
	double bin_total_xray_avgE[max_bins];
	double bin_total_xray_Tcont1[max_bins];
	double bin_total_xray_Tcont2[max_bins];
	double total_xray_Fcont;
	double total_xray_Fline;
	double total_xray_avgE;
	double total_xray_Tcont1;
	double total_xray_Tcont2;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double Tex;
	double Tex_r[num_radii];
	double bin_total_exray_Fcont[max_bins];
	double bin_total_exray_Fline[max_bins];
	double bin_total_exray_avgE[max_bins];
	double bin_total_exray_Tcont1[max_bins];
	double bin_total_exray_Tcont2[max_bins];
	double total_exray_Fcont;
	double total_exray_Fline;
	double total_exray_avgE;
	double total_exray_Tcont1;
	double total_exray_Tcont2;
#endif 

	FILE *rlist[num_radii];
	FILE *input, *stellar_input;

	char filename[256];
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char prefix[128];
	char suffix[128];
	FILE *blist;
	FILE *bszlist;
	FILE *btxlist;
	FILE *bmpro;
	FILE *bvpro;
	FILE *bzpro;
	FILE *bgpro;
#ifdef MASS_WEIGHTED
	FILE *bgpro_massweighted;
#endif /* MASS_WEIGHTED */
#ifdef COOLING
	FILE *bcpro;
#endif /* COOLING */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double Yesz, Yesz_total;
	double Yesz_r[num_radii];
	double bin_total_electron_sz_flux[max_bins];

	FILE *beszlist;
	FILE *bepro;
#endif /* ELECTRON_ION_NONEQULIBRIUM */

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

	cart_debug("Tfact = %e", Tfact );
	cart_debug("szfact = %e", szfact );

	rfact = 1000.0 * r0 * aexpn / hubble; /* code units -> proper kpc */
	vfact = v0 / aexpn;

#ifdef STARFORM
	tnewstar = 1e9 * new_star_age / t0 / (aexpn*aexpn);
#endif

#ifdef COOLING
	dtfact = 1e6 * t0 * aexpn*aexpn;
	dEfact = t0 * aexpn*aexpn/Pfact;

	fact_nH = log10( 1.12e-5 * hubble * hubble * Omega0 * ( 1.0 - Y_p ) / (aexpn*aexpn*aexpn) );
#endif

	/* set up binning */
	rlmin = log10(rbinmin/rfact);
	rlmax = log10(rbinmax/rfact);
	drl = (rlmax - rlmin)/(float)(max_bins-1);
	num_bins = (rlmax - rlmin)/drl + 1;

	cart_debug("rfact = %e", rfact );
	cart_debug("rlmin = %e", rlmin );
	cart_debug("rlmax = %e", rlmax );
	cart_debug("drl = %e", drl );
	cart_debug("num_bins = %u", num_bins );

	cart_assert( num_bins <= max_bins );

	rl[0] = 0.0;
	rr[0] = rbinmin/rfact;
	rmid[0] = 0.5*rbinmin/rfact;
	bin_volume[0] = 4.0*M_PI/3.0 * rr[0]*rr[0]*rr[0];
	bin_volume_cumulative[0] = bin_volume[0];

	for ( i = 1; i < num_bins; i++ ) {
		rl[i] = rr[i-1];
		rr[i] = pow( 10.0, rlmin + (float)i*drl );
		rmid[i] = pow( 10.0, rlmin + (float)(i-0.5)*drl );
		bin_volume[i] = 4.0*M_PI/3.0 * ( rr[i]*rr[i]*rr[i] - rl[i]*rl[i]*rl[i] );
		bin_volume_cumulative[i] = 4.0*M_PI/3.0*rr[i]*rr[i]*rr[i];
	}

	/* compute virial overdensity using Bryan & Norman (1998) */
	a3 = aexpn*aexpn*aexpn;
	E2 = Omega0 / a3 + OmegaL0 + (1.0-Omega0-OmegaL0)/(aexpn*aexpn);
	omega = Omega0 / a3 / E2;
	xz = omega - 1.0;
	virial_overdensity = ( 18.0*M_PI*M_PI + 82.0*xz - 39.0*xz*xz ) / ( 1 + xz );
	dlout = log10(virial_overdensity);

	cart_debug("E2 = %e", E2 );
	cart_debug("omega = %e", omega );
	cart_debug("xz = %e", xz );

	sprintf( suffix, "_a%6.4f.dat", aexpn );
	sprintf( prefix, "%s/h_", analysis_directory );

	/* convert overdensities and open files */
	for ( i = 0; i < num_radii; i++ ) {
		if ( radii_units[i] == 'c' ) {
			radii_overdensity[i] = radii_delta[i] / omega;
		} else {
			radii_overdensity[i] = radii_delta[i];
		}

		sprintf( filename, "%sblist_%s%s", prefix, radii_label[i], suffix );
		rlist[i] = fopen( filename, "w" );
		if ( rlist[i] == NULL ) {
			cart_error( "Unable to open %s for writing", filename );
		}
	}

	/* open output files */
	sprintf( filename, "%sblist%s", prefix, suffix );
	blist = fopen( filename, "w" );
	if ( blist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	sprintf( filename, "%sbmpro%s", prefix, suffix );
	bmpro = fopen( filename, "w" );
	if ( bmpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

#ifdef HYDRO
	sprintf( filename, "%sszlist%s", prefix, suffix );
	bszlist = fopen( filename, "w" );
	if ( bszlist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

#ifdef ANALYSIS_XRAY
	sprintf( filename, "%stxlist%s", prefix, suffix );
	btxlist = fopen( filename, "w" );
	if ( btxlist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* ANALYSIS_XRAY */

	sprintf( filename, "%sbgpro%s", prefix, suffix );
	bgpro = fopen( filename, "w" );
	if ( bgpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

#ifdef MASS_WEIGHTED
	sprintf( filename, "%sbgpro_massweighted%s", prefix, suffix );
	bgpro_massweighted = fopen( filename, "w" );
	if ( bgpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* MASS_WEIGHTED */

#ifdef COOLING
	sprintf( filename, "%sbcpro%s", prefix, suffix );
	bcpro = fopen( filename, "w" );
	if ( bgpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* COOLING */

#ifdef STARFORM
	sprintf( filename, "%sbzpro%s", prefix, suffix );
	bzpro = fopen( filename, "w" );
	if ( bzpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}
#endif /* STARFORM */
#endif /* HYDRO */

	sprintf( filename, "%sbvpro%s", prefix, suffix );
	bvpro = fopen( filename, "w" );
	if ( bvpro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	/* put headers on files */
	fprintf( blist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( blist, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( blist, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
	fprintf( blist, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( blist, "# Columns:\n");
	fprintf( blist, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( blist, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( blist, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( blist, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

#ifdef HYDRO
	fprintf( bszlist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bszlist, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bszlist, "# Overdensities: vir (%.3f)", virial_overdensity );

	for ( i = 0; i < num_radii; i++ ) {
		fprintf( bszlist, ", %s", radii_label[i] );
	}

	fprintf( bszlist, "\n");

	fprintf( bszlist, "# Columns:\n" );
	fprintf( bszlist, "# id Mg Mdm Mtotal Ysz [Mpc^2] (for each overdensity)\n" );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	sprintf( filename, "%selectron_szlist%s", prefix, suffix );
	beszlist = fopen( filename, "w" );
	if ( beszlist == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	fprintf( beszlist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( beszlist, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( beszlist, "# Overdensities: vir (%.3f)", virial_overdensity );

	for ( i = 0; i < num_radii; i++ ) {
		fprintf( beszlist, ", %s", radii_label[i] );
	}

	fprintf( beszlist, "\n");

	fprintf( beszlist, "# Columns:\n" );
	fprintf( beszlist, "# id Ysz Yesz [Mpc^2] (for each overdensity)\n" );
#endif

#ifdef ANALYSIS_XRAY
	fprintf( btxlist, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( btxlist, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( btxlist, "# Overdensities: vir (%.3f)", virial_overdensity );

	for ( i = 0; i < num_radii; i++ ) {
		fprintf( btxlist, ", %s", radii_label[i] );
	}

	fprintf( btxlist, "\n");

	fprintf( btxlist, "# Columns:\n" );
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	fprintf( btxlist, "# id Tx Tex (for each overdensity)\n" );	
#else
	fprintf( btxlist, "# id Mg Mdm Mtotal Tx (for each overdensity)\n" );
#endif
#endif /* ANALYSIS_XRAY */
#endif /* HYDRO */

	/* density profiles */
	fprintf( bmpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bmpro, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bmpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
	fprintf( bmpro, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( bmpro, "# Header Columns:\n");
	fprintf( bmpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( bmpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( bmpro, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( bmpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

	fprintf( bmpro, "# Profile Columns:\n");
	fprintf( bmpro, "# rmid rr [/h kpc] Mdm Mg Mcg M* M*new (< rr) [/h Msolar]\n");
	fprintf( bmpro, "#########################################################################\n");

	/* velocity profiles */
	fprintf( bvpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bvpro, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bvpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
	fprintf( bvpro, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( bvpro, "# Header Columns:\n");
	fprintf( bvpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( bvpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( bvpro, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( bvpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

	fprintf( bvpro, "# Profile Columns:\n");
	fprintf( bvpro, "# rmid rr [/h kpc] Vrms_dm Vrms_gas Vrms_* Vrms_*new [km/s]\n");
	fprintf( bvpro, "# Vc_dm Vc_gas Vc_* Vc_total\n");
	fprintf( bvpro, "#########################################################################\n");

#ifdef HYDRO
	fprintf( bgpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bgpro, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bgpro, "# Tcold = %eK, new_star_age = %.2f Gyr\n", Tcold, new_star_age );
	fprintf( bgpro, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( bgpro, "# Header Columns:\n");
	fprintf( bgpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( bgpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( bgpro, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( bgpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

	fprintf( bgpro, "# Profile Columns:\n");
	fprintf( bgpro, "# rmid rr [/h kpc] Mg Mcg [/h Msolar] T [K] P [ergs cm^{-3}] S [keV cm^2]\n");
	fprintf( bgpro, "#########################################################################\n");

#ifdef MASS_WEIGHTED 
    fprintf( bgpro_massweighted, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
    fprintf( bgpro_massweighted, "# Monte carlo points per cell = %u\n", points_per_cell );
    fprintf( bgpro_massweighted, "# Tcold = %eK, new_star_age = %.2f Gyr\n", Tcold, new_star_age );
    fprintf( bgpro_massweighted, "# virial overdensity = %.3f\n", virial_overdensity );

    fprintf( bgpro_massweighted, "# Header Columns:\n");
    fprintf( bgpro_massweighted, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
    fprintf( bgpro_massweighted, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
    fprintf( bgpro_massweighted, "# vmax [km/s] rmax [/h kpc]\n" );
    fprintf( bgpro_massweighted, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

    fprintf( bgpro_massweighted, "# Profile Columns:\n");
    fprintf( bgpro_massweighted, "# rmid rr [/h kpc] Mg Mcg [/h Msolar] T [K] P [ergs cm^{-3}] S [keV cm^-2]\n");
	fprintf( bgpro_massweighted, "#########################################################################\n");
#endif /* MASS_WEIGHTED */

#ifdef COOLING
	fprintf( bcpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bcpro, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bcpro, "# Tcold = %eK, new_star_age = %.2f Gyr\n", Tcold, new_star_age );
	fprintf( bcpro, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( bcpro, "# Header Columns:\n");
	fprintf( bcpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( bcpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( bcpro, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( bcpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

	fprintf( bcpro, "# Profile Columns:\n");
	fprintf( bcpro, "# rmid rr [/h kpc] Mg Mcg [/h Msolar] dE_cool [ergs s^-1] <tcool> [Myr] <Zg_II> <Zg_Ia>\n");
	fprintf( bcpro, "#########################################################################\n");
#endif /* COOLING */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	sprintf( filename, "%sbepro%s", prefix, suffix );
	bepro = fopen( filename, "w" );
	if ( bepro == NULL ) {
		cart_error("Unable to open %s for writing.", filename );
	}

	fprintf( bepro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bepro, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bepro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
	fprintf( bepro, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( bepro, "# Header Columns:\n");
	fprintf( bepro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( bepro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( bepro, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( bepro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

	fprintf( bepro, "# Profile Columns:\n");
	fprintf( bepro, "# rmid rr [/h kpc] Mg [/h Msolar] T [K] Te [K] Ysz Ysz_electron [Mpc^2] n_e t_ei Pgas Pe\n");
	fprintf( bepro, "#########################################################################\n");
#endif

#ifdef STARFORM
	/* star and metallicity profiles */
	fprintf( bzpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( bzpro, "# Monte carlo points per cell = %u\n", points_per_cell );
	fprintf( bzpro, "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
	fprintf( bzpro, "# virial overdensity = %.3f\n", virial_overdensity );

	fprintf( bzpro, "# Header Columns:\n");
	fprintf( bzpro, "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
	fprintf( bzpro, "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rvir] Mvir (from halo finder)\n");
	fprintf( bzpro, "# vmax [km/s] rmax [/h kpc]\n" );
	fprintf( bzpro, "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

	fprintf( bzpro, "# Profile Columns:\n");
	fprintf( bzpro, "# rmid rr [/h kpc] Zg_II Zg_Ia Z*_II Z*_Ia Z*new_II Z*new_Ia <t*> [Gyr]\n");
	fprintf( bzpro, "#########################################################################\n");
#endif /* STARFORM */
#endif /* HYDRO */

	for ( i = 0; i < num_radii; i++ ) {
		fprintf( rlist[i], "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
		fprintf( rlist[i], "# Monte carlo points per cell = %u\n", points_per_cell );
		fprintf( rlist[i], "# Tcold = %fK, new_star_age = %f Gyr\n", Tcold, new_star_age );
		fprintf( rlist[i], "# overdensity = %s, %.3f\n", radii_label[i], radii_delta[i] );

		if ( radii_units[i] == 'r' ) {
			fprintf( rlist[i], "# Quantities at fixed radius rout\n" );
		} else if ( radii_units[i] == 'c' ) {
			fprintf( rlist[i], "# Quantities at fixed overdensity %f with respect to critical\n", radii_delta[i] );
		} else if ( radii_units[i] == 'm' ) {
			fprintf( rlist[i], "# Quantities at fixed matter overdensity %f\n", radii_delta[i] );
		}

		fprintf( rlist[i], "# Columns:\n");
		fprintf( rlist[i], "# id min(rout,rmax) min(rvir,rmax) [/h kpc]\n" );
		fprintf( rlist[i], "# Mg Mcg M* M*new Mb Mdm Mtotal [/h Msolar,< rout] Mvir [/h Msolar]\n");
		fprintf( rlist[i], "# vmax [km/s] rmax [/h kpc]\n" );
		fprintf( rlist[i], "# <Zg_II> <Zg_Ia> <Z*_II> <Z*_Ia> <Z*new_II> <Z*new_Ia> <t*> [Gyr]\n");

		fflush(rlist[i]);
	}

	/* process all halos */
	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		profile = cart_alloc( sizeof(profile_struct) );

		/* clear out bins */
		for ( bin = 0; bin < num_bins; bin++ ) {
#ifdef PARTICLES
			profile->bin_dark_num[bin] = 0;
			profile->bin_dark_mass[bin] = 0.0;
			profile->bin_dark_momentum[0][bin] = 0.0;
			profile->bin_dark_momentum[1][bin] = 0.0;
			profile->bin_dark_momentum[2][bin] = 0.0;
			profile->bin_dark_vrms[bin] = 0.0;
#ifdef STARFORM 
			profile->bin_star_num[bin] = 0;
			profile->bin_star_mass[bin] = 0.0;
			profile->bin_star_momentum[0][bin] = 0.0;
			profile->bin_star_momentum[1][bin] = 0.0;
			profile->bin_star_momentum[2][bin] = 0.0;
			profile->bin_star_vrms[bin] = 0.0;
			profile->bin_star_age[bin] = 0.0;
			profile->bin_star_metallicity_II[bin] = 0.0;
			profile->bin_star_metallicity_Ia[bin] = 0.0;
			profile->bin_new_star_num[bin] = 0;
			profile->bin_new_star_mass[bin] = 0.0;
			profile->bin_new_star_momentum[0][bin] = 0.0;
			profile->bin_new_star_momentum[1][bin] = 0.0;
			profile->bin_new_star_momentum[2][bin] = 0.0;
			profile->bin_new_star_vrms[bin] = 0.0;
			profile->bin_new_star_metallicity_II[bin] = 0.0;
			profile->bin_new_star_metallicity_Ia[bin] = 0.0;
#endif /* STARFORM */
#endif /* PARTICLES */

			profile->bin_gas_mass[bin] = 0.0;
			profile->bin_excluded_gas_mass[bin] = 0.0;
			profile->bin_volume[bin] = 0.0;
			profile->bin_excluded_volume[bin] = 0.0;

			profile->bin_gas_velocity[0][bin] = 0.0;
			profile->bin_gas_velocity[1][bin] = 0.0;
			profile->bin_gas_velocity[2][bin] = 0.0;
			profile->bin_gas_vrms[bin] = 0.0;
			profile->bin_cold_gas_mass[bin] = 0.0;

			profile->bin_gas_temperature[bin] = 0.0;
			profile->bin_gas_entropy[bin] = 0.0;
			profile->bin_gas_pressure[bin] = 0.0;

#ifdef MASS_WEIGHTED
			profile->bin_massweighted_temperature[bin] = 0.0;
			profile->bin_massweighted_entropy[bin] = 0.0;
			profile->bin_massweighted_pressure[bin] = 0.0;
#endif /* MASS_WEIGHTED */

			profile->bin_gas_coolingrate[bin] = 0.0;
			profile->bin_gas_tcool[bin] = 0.0;
			profile->bin_sz_flux[bin] = 0.0;

			profile->bin_gas_metallicity_II[bin] = 0.0;
			profile->bin_gas_metallicity_Ia[bin] = 0.0;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			profile->bin_electron_pressure[bin] = 0.0;
			profile->bin_electron_temperature[bin] = 0.0;
			profile->bin_electron_sz_flux[bin] = 0.0;
			profile->bin_electron_density[bin] = 0.0;
			profile->bin_electron_tei[bin] = 0.0;
#endif

#ifdef ANALYSIS_XRAY
			profile->bin_xray_Fcont[bin] = 0.0;
			profile->bin_xray_Fline[bin] = 0.0;
			profile->bin_xray_avgE[bin] = 0.0;
			profile->bin_xray_Tcont1[bin] = 0.0;
			profile->bin_xray_Tcont2[bin] = 0.0;

#ifdef ELECTRON_ION_NONEQUILIBRIUM 
			profile->bin_exray_Fcont[bin] = 0.0;
			profile->bin_exray_Fline[bin] = 0.0;
			profile->bin_exray_avgE[bin] = 0.0;
			profile->bin_exray_Tcont1[bin] = 0.0;
			profile->bin_exray_Tcont2[bin] = 0.0;
#endif
#endif /* ANALYSIS_XRAY */
		}

		halos->list[ihalo].analysis_radius = 1.05*rbinmax/(1000.0);

		profile->rlmin = rlmin;
		profile->drl = drl;
		profile->num_bins = num_bins;

		halos->list[ihalo].data = profile;
	}

#ifdef PARTICLES
	/* open particle file */
	sprintf( filename1, "%s/PMcrda%06.4f.DAT", output_directory, aexpn );
	sprintf( filename2, "%s/PMcrs_indexed_a%06.4f.DAT", output_directory, aexpn );
	sprintf( filename3, "%s/stars_indexed_a%06.4f.dat", output_directory, aexpn );

	read_indexed_particles( filename1, filename2, filename3, halos, subhalos, particle_callback );

	cart_debug("done assigning particles");
#endif /* PARTICLES */

#ifdef HYDRO
	sprintf(filename, "%s/%s_a%6.4f.grid", output_directory, jobname, aexpn );
	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s for reading!\n", filename );
	}

	read_indexed_grid( filename, halos, subhalos, cell_callback );

	cart_debug("done assigning gas to grid"); fflush(stdout);
#endif /* HYDRO */

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		profile = (profile_struct *)halos->list[ihalo].data;

#ifdef HYDRO
		/* process bins */
		for ( bin = 0; bin < num_bins; bin++ ) {
			if ( profile->bin_gas_mass[bin] > 0.0 ) {
				inverse_mass = 1.0 / profile->bin_excluded_gas_mass[bin];

#ifdef MASS_WEIGHTED
				profile->bin_massweighted_pressure[bin] *= inverse_mass;
				profile->bin_massweighted_temperature[bin] *= inverse_mass;
				profile->bin_massweighted_entropy[bin] *= inverse_mass;
#endif /* MASS_WEIGHTED */

				profile->bin_gas_pressure[bin] /= profile->bin_excluded_volume[bin];
				profile->bin_gas_temperature[bin] /= profile->bin_excluded_volume[bin];
				profile->bin_gas_entropy[bin] /= profile->bin_excluded_volume[bin];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
				profile->bin_electron_pressure[bin] /= profile->bin_excluded_volume[bin];
				profile->bin_electron_temperature[bin] /= profile->bin_excluded_volume[bin];
				profile->bin_electron_density[bin] /= profile->bin_excluded_volume[bin];
				profile->bin_electron_tei[bin] /= profile->bin_excluded_volume[bin];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

				profile->bin_gas_tcool[bin] *= inverse_mass;
				profile->bin_gas_coolingrate[bin] /= profile->bin_excluded_volume[bin];

				profile->bin_gas_metallicity_II[bin] *= inverse_mass;
				profile->bin_gas_metallicity_Ia[bin] *= inverse_mass;

				profile->bin_gas_velocity[0][bin] *= inverse_mass;
				profile->bin_gas_velocity[1][bin] *= inverse_mass;
				profile->bin_gas_velocity[2][bin] *= inverse_mass;
				vmean = profile->bin_gas_velocity[0][bin] *
					profile->bin_gas_velocity[0][bin] +
					profile->bin_gas_velocity[1][bin] *
					profile->bin_gas_velocity[1][bin] +
					profile->bin_gas_velocity[2][bin] *
					profile->bin_gas_velocity[2][bin];
				profile->bin_gas_vrms[bin] = 
					sqrt(fabs(profile->bin_gas_vrms[bin]*inverse_mass - vmean))*vfact;
				profile->bin_gas_velocity[0][bin] *= vfact;
				profile->bin_gas_velocity[1][bin] *= vfact;
				profile->bin_gas_velocity[2][bin] *= vfact;
			}
		}	
#endif /* HYDRO */

		/* compute cumulative quantities */
		total_dark_mass = 0.0;
		total_star_mass = 0.0;
		total_new_star_mass = 0.0;
		total_gas_mass = 0.0;
		total_cold_gas_mass = 0.0;

		Ysz_total = 0.0;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		Yesz_total = 0.0;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

		avg_gas_metallicity_II = 0.0;
		avg_gas_metallicity_Ia = 0.0;
		avg_star_metallicity_II = 0.0;
		avg_star_metallicity_Ia = 0.0;
		avg_new_star_metallicity_II = 0.0;
		avg_new_star_metallicity_Ia = 0.0;
		avg_star_age = 0.0;

		avg_new_star_mass = 0.0;
		avg_star_mass = 0.0;
		avg_gas_mass = 0.0;

#ifdef ANALYSIS_XRAY
		total_xray_Fcont = 0.0;
		total_xray_Fline = 0.0;
		total_xray_avgE = 0.0;
		total_xray_Tcont1 = 0.0;
		total_xray_Tcont2 = 0.0;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
		total_exray_Fcont = 0.0;
		total_exray_Fline = 0.0;
		total_exray_avgE = 0.0;
		total_exray_Tcont1 = 0.0;
		total_exray_Tcont2 = 0.0;
#endif
#endif

		irvir = 0;
		irflag = 0;

		rvir = 0.0;
		rout = 0.0;

		for ( i = 0; i < num_radii; i++ ) {
			avg_gas_metallicity_II_r[i] = 0.0;
			avg_gas_metallicity_Ia_r[i] = 0.0;
			avg_star_metallicity_II_r[i] = 0.0;
			avg_star_metallicity_Ia_r[i] = 0.0;
			avg_new_star_metallicity_II_r[i] = 0.0;
			avg_new_star_metallicity_Ia_r[i] = 0.0;
			avg_star_age_r[i] = 0.0;

			avg_new_star_mass_r[i] = 0.0;
			avg_star_mass_r[i] = 0.0;
			avg_gas_mass_r[i] = 0.0;

			iflag_r[i] = 0;
			delta_r[i] = 0.0;
		}

		vmax = 0.0;

		for ( bin = 0; bin < num_bins; bin++ ) {
			total_dark_mass += profile->bin_dark_mass[bin];
			total_star_mass += profile->bin_star_mass[bin];
			total_new_star_mass += profile->bin_new_star_mass[bin];
			total_gas_mass += profile->bin_gas_mass[bin];
			total_cold_gas_mass += profile->bin_cold_gas_mass[bin];

			Ysz_total += profile->bin_sz_flux[bin];

			bin_total_dark_mass[bin] = total_dark_mass;
			bin_total_star_mass[bin] = total_star_mass;
			bin_total_new_star_mass[bin] = total_new_star_mass;
			bin_total_gas_mass[bin] = total_gas_mass;
			bin_total_cold_gas_mass[bin] = total_cold_gas_mass;
			bin_total_sz_flux[bin] = Ysz_total*rfact*rfact/1e6;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			Yesz_total += profile->bin_electron_sz_flux[bin];
			bin_total_electron_sz_flux[bin] = Yesz_total*rfact*rfact/1e6;
#endif

#ifdef ANALYSIS_XRAY
			total_xray_Fcont += profile->bin_xray_Fcont[bin];
			total_xray_Fline += profile->bin_xray_Fline[bin];
			total_xray_avgE += profile->bin_xray_avgE[bin];
			total_xray_Tcont1 += profile->bin_xray_Tcont1[bin];
			total_xray_Tcont2 += profile->bin_xray_Tcont2[bin];

			bin_total_xray_Fcont[bin] = total_xray_Fcont;
			bin_total_xray_Fline[bin] = total_xray_Fline;
			bin_total_xray_avgE[bin] = total_xray_avgE;
			bin_total_xray_Tcont1[bin] = total_xray_Tcont1;
			bin_total_xray_Tcont2[bin] = total_xray_Tcont2;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			total_exray_Fcont += profile->bin_exray_Fcont[bin];
			total_exray_Fline += profile->bin_exray_Fline[bin];
			total_exray_avgE += profile->bin_exray_avgE[bin];
			total_exray_Tcont1 += profile->bin_exray_Tcont1[bin];
			total_exray_Tcont2 += profile->bin_exray_Tcont2[bin];

			bin_total_exray_Fcont[bin] = total_exray_Fcont; 
			bin_total_exray_Fline[bin] = total_exray_Fline;
			bin_total_exray_avgE[bin] = total_exray_avgE;
			bin_total_exray_Tcont1[bin] = total_exray_Tcont1;
			bin_total_exray_Tcont2[bin] = total_exray_Tcont2;
#endif
#endif /* ANALYSIS_XRAY */

			if ( irflag == 0 ) {
				avg_gas_metallicity_II += profile->bin_gas_metallicity_II[bin] * 
					profile->bin_gas_mass[bin];
				avg_gas_metallicity_Ia += profile->bin_gas_metallicity_Ia[bin] * 
					profile->bin_gas_mass[bin];
				avg_star_metallicity_II += profile->bin_star_metallicity_II[bin];
				avg_star_metallicity_Ia += profile->bin_star_metallicity_Ia[bin];
				avg_new_star_metallicity_II += profile->bin_new_star_metallicity_II[bin];
				avg_new_star_metallicity_Ia += profile->bin_new_star_metallicity_Ia[bin];
				avg_star_age += profile->bin_star_age[bin];

				avg_gas_mass += profile->bin_gas_mass[bin];
				avg_new_star_mass += profile->bin_new_star_mass[bin];
				avg_star_mass += profile->bin_star_mass[bin];
			}

			/* find maximum circular velocity within min(rvir,rhalo) */
			if ( bin > 0 && irflag == 0 && rr[bin] < halos->list[ihalo].rhalo) {
				if ( sqrt( (	bin_total_dark_mass[bin] +
								bin_total_gas_mass[bin] +
								bin_total_star_mass[bin] ) / rr[bin] ) >= vmax ) {
					rmax = rr[bin];
					vmax = sqrt( ( 	bin_total_dark_mass[bin] +
								bin_total_gas_mass[bin] +
								bin_total_star_mass[bin] ) / rr[bin] );
				}
			}

			if ( bin == 0 ) {
				dbi1 = 0.0;
			} else {
				dbi1 = ( bin_total_dark_mass[bin-1] + bin_total_star_mass[bin-1] + 
						bin_total_gas_mass[bin-1] ) / bin_volume_cumulative[bin-1];
			}

			dbi2 = (bin_total_dark_mass[bin] + bin_total_star_mass[bin] + 
					bin_total_gas_mass[bin]) / bin_volume_cumulative[bin];

			/* compute virial radius */
			if ( bin > 0 && ( dbi1 >= virial_overdensity && dbi2 < virial_overdensity ) ) {
				rrl = log10(rr[bin]);
				rll = log10(rl[bin]);
				rri = 1.0/(rrl-rll);
				dlbi1 = log10(dbi1);
				dlbi2 = log10(dbi2);
				irvir = 1;
				rvir = pow( 10.0, (dlout*(rrl-rll) + rll*dlbi2 - 
							rrl*dlbi1)/(dlbi2-dlbi1));
			}

			if ( bin > 0 && irflag == 0 &&
					( ( rl[bin] <= rmass && rr[bin] > rmass ) || 
					  ( dbi1 >= virial_overdensity && dbi2 < virial_overdensity ) ) ) {

				irflag = 1;
				rrl = log10(rr[bin]);
				rll = log10(rl[bin]);
				rri = 1.0/(rrl-rll);

				if ( dbi1 >= virial_overdensity && dbi2 < virial_overdensity ) {
					dlbi1 = log10(dbi1);
					dlbi2 = log10(dbi2);
					rout = pow( 10.0, (dlout*(rrl-rll) + rll*dlbi2 
								- rrl*dlbi1)/(dlbi2-dlbi1));
				} else {
					rout = rmass;
				}

				rlout = log10(rout);

				/* interpolate total quantities to rout */
				aM_gas = log_interpolate( bin_total_gas_mass, bin, rlout, rri, rll );
				aM_cold_gas = log_interpolate( bin_total_cold_gas_mass, bin, rlout, rri, rll );

#ifdef STARFORM
				aM_stars = log_interpolate( bin_total_star_mass, bin, rlout, rri, rll );
				aM_new_stars = log_interpolate( bin_total_new_star_mass, bin, rlout, rri, rll );
#else
				aM_stars = 0.0;
				aM_new_stars = 0.0;
#endif /* STARFORM */

				aM_dark = log_interpolate( bin_total_dark_mass, bin, rlout, rri, rll );

				Ysz = log_interpolate( bin_total_sz_flux, bin, rlout, rri, rll );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
				Yesz = log_interpolate( bin_total_electron_sz_flux, bin, rlout, rri, rll );
#endif

#ifdef ANALYSIS_XRAY
				Fcont  = log_interpolate( bin_total_xray_Fcont,  bin, rlout, rri, rll );
				Fline  = log_interpolate( bin_total_xray_Fline,  bin, rlout, rri, rll );
				avgE   = log_interpolate( bin_total_xray_avgE,   bin, rlout, rri, rll );
				Tcont1 = log_interpolate( bin_total_xray_Tcont1, bin, rlout, rri, rll );
				Tcont2 = log_interpolate( bin_total_xray_Tcont2, bin, rlout, rri, rll );

				/* subtract off inner 0.15 R_x */
				/*
				j = 0;
				while ( j < bin ) {
					if ( rr[j] >= 0.15*rout ) {
						break;
					}
					j++;
				}

				rlout = log10(0.15*rout);
				rrl = log10(rr[j]);
				rll = log10(rl[j]);
				rri = 1.0/(rrl-rll);

				Fcont -= log_interpolate( bin_total_xray_Fcont, j, rlout, rri, rll );
				Fline -= log_interpolate( bin_total_xray_Fline, j, rlout, rri, rll );
				avgE -= log_interpolate( bin_total_xray_avgE, j, rlout, rri, rll );
				Tcont1 -= log_interpolate( bin_total_xray_Tcont1, j, rlout, rri, rll );
				Tcont2 -= log_interpolate( bin_total_xray_Tcont2, j, rlout, rri, rll );
				*/

				Tcont = Tcont1/Tcont2;
				avgE /= Fline;
				Tline = xray_calibrated_line_temperature(avgE);

				f_line = Fline / ( Fline + Fcont );
				xx = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
					exp( -pow(f_line/delta_xray_2,8.0) );

				Tx = xx*Tcont + (1.0-xx)*Tline;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
				Fcont  = log_interpolate( bin_total_exray_Fcont,  bin, rlout, rri, rll );
				Fline  = log_interpolate( bin_total_exray_Fline,  bin, rlout, rri, rll );
				avgE   = log_interpolate( bin_total_exray_avgE,   bin, rlout, rri, rll );
				Tcont1 = log_interpolate( bin_total_exray_Tcont1, bin, rlout, rri, rll );
				Tcont2 = log_interpolate( bin_total_exray_Tcont2, bin, rlout, rri, rll );

				/* subtract off inner 0.15 R_x 
				   j = 0;
				   while ( j < bin ) {
				   if ( rr[j] >= 0.15*rout ) {
				   break;
				   }
				   j++;
				   }

				   rlout = log10(0.15*rout);
				   rrl = log10(rr[j]);
				   rll = log10(rl[j]);
				   rri = 1.0/(rrl-rll);

				   Fcont -= log_interpolate( bin_total_exray_Fcont, j, rlout, rri, rll );
				   Fline -= log_interpolate( bin_total_exray_Fline, j, rlout, rri, rll );
				   avgE -= log_interpolate( bin_total_exray_avgE, j, rlout, rri, rll );
				   Tcont1 -= log_interpolate( bin_total_exray_Tcont1, j, rlout, rri, rll );
				   Tcont2 -= log_interpolate( bin_total_exray_Tcont2, j, rlout, rri, rll );
				 */

				Tcont = Tcont1/Tcont2;
				avgE /= Fline;
				Tline = xray_calibrated_line_temperature(avgE);

				f_line = Fline / ( Fline + Fcont );
				xx = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
					exp( -pow(f_line/delta_xray_2,8.0) );

				Tex = xx*Tcont + (1.0-xx)*Tline;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */

				/* convert to output units */
				aM_gas *= aM0 * hubble;
				aM_cold_gas *= aM0 * hubble;
				aM_stars *= aM0 * hubble;
				aM_new_stars *= aM0 * hubble;
				aM_dark *= aM0 * hubble;
			}

			for ( i = 0; i < num_radii; i++ ) {
				if ( iflag_r[i] == 0 ) {
					avg_gas_metallicity_II_r[i] += profile->bin_gas_metallicity_II[bin] * 
						profile->bin_gas_mass[bin];
					avg_gas_metallicity_Ia_r[i] += profile->bin_gas_metallicity_Ia[bin] * 
						profile->bin_gas_mass[bin];
					avg_star_metallicity_II_r[i] += profile->bin_star_metallicity_II[bin];
					avg_star_metallicity_Ia_r[i] += profile->bin_star_metallicity_Ia[bin];
					avg_new_star_metallicity_II_r[i] += profile->bin_new_star_metallicity_II[bin];
					avg_new_star_metallicity_Ia_r[i] += profile->bin_new_star_metallicity_Ia[bin];

					avg_star_age_r[i] += profile->bin_star_age[bin];
					avg_gas_mass_r[i] += profile->bin_gas_mass[bin];
					avg_new_star_mass_r[i] += profile->bin_new_star_mass[bin];
					avg_star_mass_r[i] += profile->bin_star_mass[bin];
				}

				if ( bin > 0 && iflag_r[i] == 0 &&
						( ( i > 0 && dbi1 >= radii_overdensity[i] && 
							dbi2 < radii_overdensity[i] ) ||
						  ( i == 0 && rl[bin] < halos->list[ihalo].rhalo && 
							rr[bin] >= halos->list[ihalo].rhalo ) ) ) {

					iflag_r[i] = 1;
					rrl = log10(rr[bin]);
					rll = log10(rl[bin]);
					rri = 1.0/(rrl-rll);

					if ( i == 0 ) {
						delta_r[i] = halos->list[ihalo].rhalo;
					} else {
						dlbi1 = log10(dbi1);
						dlbi2 = log10(dbi2);
						delta_r[i] = pow( 10.0, (log10(radii_overdensity[i])*(rrl-rll) + rll*dlbi2
									- rrl*dlbi1)/(dlbi2-dlbi1));
					}

					rlout = log10(delta_r[i]);

					/* interpolate total quantities to rlout */
					aM_gas_r[i] = log_interpolate( bin_total_gas_mass, bin, rlout, rri, rll );
					aM_cold_gas_r[i] = log_interpolate( bin_total_cold_gas_mass, bin, rlout, rri, rll );

#ifdef STARFORM
					aM_stars_r[i] = log_interpolate( bin_total_star_mass, bin, rlout, rri, rll );
					aM_new_stars_r[i] = log_interpolate( bin_total_new_star_mass, bin, rlout, rri, rll );
#else
					aM_stars_r[i] = 0.0;
					aM_new_stars_r[i] = 0.0;
#endif /* STARFORM */

					aM_dark_r[i] = log_interpolate( bin_total_dark_mass, bin, rlout, rri, rll );

					Ysz_r[i] = log_interpolate( bin_total_sz_flux, bin, rlout, rri, rll );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					Yesz_r[i] = log_interpolate( bin_total_electron_sz_flux, bin, rlout, rri, rll );
#endif

#ifdef ANALYSIS_XRAY
					Fcont = log_interpolate( bin_total_xray_Fcont, bin, rlout, rri, rll );
					Fline = log_interpolate( bin_total_xray_Fline, bin, rlout, rri, rll );
					avgE = log_interpolate( bin_total_xray_avgE, bin, rlout, rri, rll );
					Tcont1 = log_interpolate( bin_total_xray_Tcont1, bin, rlout, rri, rll );
					Tcont2 = log_interpolate( bin_total_xray_Tcont2, bin, rlout, rri, rll );

					/* subtract off inner 0.1 R_x 
					rlout = 0.15*delta_r[i];
					j = 0;
					while ( j < bin ) {
						if ( rr[j] >= rlout ) {
							break;
						}
						j++;
					}

					rlout = log10(rlout);
					rrl = log10(rr[j]);
					rll = log10(rl[j]);
					rri = 1.0/(rrl-rll);

					Fcont -= log_interpolate( bin_total_xray_Fcont, j, rlout, rri, rll );
					Fline -= log_interpolate( bin_total_xray_Fline, j, rlout, rri, rll );
					avgE -= log_interpolate( bin_total_xray_avgE, j, rlout, rri, rll );
					Tcont1 -= log_interpolate( bin_total_xray_Tcont1, j, rlout, rri, rll );
					Tcont2 -= log_interpolate( bin_total_xray_Tcont2, j, rlout, rri, rll );
					*/

					Tcont = Tcont1/Tcont2;
					avgE /= Fline;
					Tline = xray_calibrated_line_temperature(avgE);

					f_line = Fline / ( Fline + Fcont );
					xx = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
						exp( -pow(f_line/delta_xray_2,8.0) );

					Tx_r[i] = xx*Tcont + (1.0-xx)*Tline;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					Fcont = log_interpolate( bin_total_exray_Fcont, bin, rlout, rri, rll );
					Fline = log_interpolate( bin_total_exray_Fline, bin, rlout, rri, rll );
					avgE = log_interpolate( bin_total_exray_avgE, bin, rlout, rri, rll );
					Tcont1 = log_interpolate( bin_total_exray_Tcont1, bin, rlout, rri, rll );
					Tcont2 = log_interpolate( bin_total_exray_Tcont2, bin, rlout, rri, rll );

					/* subtract off inner 0.1 R_x 
					   rlout = 0.15*delta_r[i];
					   j = 0;
					   while ( j < bin ) {
					   if ( rr[j] >= rlout ) {
					   break;
					   }
					   j++;
					   }

					   rlout = log10(rlout);
					   rrl = log10(rr[j]);
					   rll = log10(rl[j]);
					   rri = 1.0/(rrl-rll);

					   Fcont -= log_interpolate( bin_total_exray_Fcont, j, rlout, rri, rll );
					   Fline -= log_interpolate( bin_total_exray_Fline, j, rlout, rri, rll );
					   avgE -= log_interpolate( bin_total_exray_avgE, j, rlout, rri, rll );
					   Tcont1 -= log_interpolate( bin_total_exray_Tcont1, j, rlout, rri, rll );
					   Tcont2 -= log_interpolate( bin_total_exray_Tcont2, j, rlout, rri, rll );
					 */

					Tcont = Tcont1/Tcont2;
					avgE /= Fline;
					Tline = xray_calibrated_line_temperature(avgE);

					f_line = Fline / ( Fline + Fcont );
					xx = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
						exp( -pow(f_line/delta_xray_2,8.0) );

					Tex_r[i] = xx*Tcont + (1.0-xx)*Tline;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */

					aM_gas_r[i] *= aM0 * hubble;
					aM_cold_gas_r[i] *= aM0 * hubble;
					aM_stars_r[i] *= aM0 * hubble;
					aM_new_stars_r[i] *= aM0 * hubble;
					aM_dark_r[i] *= aM0 * hubble;
				}
			}
		}

		if ( irvir == 0 ) {
			/* never reached virial radius */
			rout = rbinmax/rfact;
			rvir = rbinmax/rfact;

			aM_gas = total_gas_mass * aM0 * hubble;
			aM_cold_gas = total_cold_gas_mass * aM0 * hubble;
			aM_stars = total_star_mass * aM0 * hubble;
			aM_new_stars = total_new_star_mass * aM0 * hubble;
			aM_dark = total_dark_mass * aM0 * hubble;

			Ysz = Ysz_total;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			Yesz = Yesz_total;
#endif

#ifdef ANALYSIS_XRAY
			Fcont = bin_total_xray_Fcont[max_bins-1];
			Fline = bin_total_xray_Fline[max_bins-1];
			avgE = bin_total_xray_avgE[max_bins-1];
			Tcont1 = bin_total_xray_Tcont1[max_bins-1];
			Tcont2 = bin_total_xray_Tcont2[max_bins-1];

			Tcont = Tcont1/Tcont2;
			avgE /= Fline;
			Tline = xray_calibrated_line_temperature(avgE);

			f_line = Fline / ( Fline + Fcont );
			xx = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
				exp( -pow(f_line/delta_xray_2,8.0) );

			Tx = xx*Tcont + (1.0-xx)*Tline;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			Fcont = bin_total_exray_Fcont[max_bins-1];
			Fline = bin_total_exray_Fline[max_bins-1];
			avgE = bin_total_exray_avgE[max_bins-1];
			Tcont1 = bin_total_exray_Tcont1[max_bins-1];
			Tcont2 = bin_total_exray_Tcont2[max_bins-1];

			Tcont = Tcont1/Tcont2;
			avgE /= Fline;
			Tline = xray_calibrated_line_temperature(avgE);

			f_line = Fline / ( Fline + Fcont );
			xx = exp( -pow(f_line/delta_xray_1,2*beta_xray) ) *
				exp( -pow(f_line/delta_xray_2,8.0) );

			Tex = xx*Tcont + (1.0-xx)*Tline;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* ANALYSIS_XRAY */
		}

		/* compute average quantities */
		if ( avg_gas_mass > 0.0 ) {
			avg_gas_metallicity_II /= avg_gas_mass*Zsolar; 
			avg_gas_metallicity_Ia /= avg_gas_mass*Zsolar;
		}

		if ( avg_star_mass > 0.0 ) {
			avg_star_metallicity_II /= avg_star_mass*Zsolar;
			avg_star_metallicity_Ia /= avg_star_mass*Zsolar;
			avg_star_age = age_t( avg_star_age / avg_star_mass );
		}

		if ( avg_new_star_mass > 0.0 ) {
			avg_new_star_metallicity_II /= avg_new_star_mass*Zsolar;
			avg_new_star_metallicity_Ia /= avg_new_star_mass*Zsolar;
		}

		aM_baryons = aM_stars + aM_gas;
		aM_total = aM_baryons + aM_dark;

		rdout = rout * rfact * hubble;
		rvdout = rvir * rfact * hubble;

		rmax *= 1000.0*r0;
		vmax *= 2.07498e-3 * sqrt(aM0/(1000.*r0/hubble));

		/* output */
		fprintf( blist, "%u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n", 
				halos->list[ihalo].id, rdout, rvdout, aM_gas, 
				aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir, vmax, rmax, 
				avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia, 
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );
		fflush(blist);

		for ( i = 0; i < num_radii; i++ ) {
			/* compute average quantities */
			if ( avg_gas_mass_r[i] > 0.0 ) {
				avg_gas_metallicity_II_r[i] /= avg_gas_mass_r[i]*Zsolar;
				avg_gas_metallicity_Ia_r[i] /= avg_gas_mass_r[i]*Zsolar;
			}

			if ( avg_star_mass_r[i] > 0.0 ) {
				avg_star_metallicity_II_r[i] /= avg_star_mass_r[i]*Zsolar;
				avg_star_metallicity_Ia_r[i] /= avg_star_mass_r[i]*Zsolar;
				avg_star_age_r[i] = age_t( avg_star_age_r[i] / avg_star_mass_r[i] );
			}

			if ( avg_new_star_mass_r[i] > 0.0 ) {
				avg_new_star_metallicity_II_r[i] /= avg_new_star_mass_r[i]*Zsolar;
				avg_new_star_metallicity_Ia_r[i] /= avg_new_star_mass_r[i]*Zsolar;
			}

			aM_baryons_r[i] = aM_stars_r[i] + aM_gas_r[i];
			aM_total_r[i] = aM_baryons_r[i] + aM_dark_r[i];

			fprintf( rlist[i], "%u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
					halos->list[ihalo].id, 
					delta_r[i]*rfact*hubble, rvdout,
					aM_gas_r[i], aM_cold_gas_r[i], aM_stars_r[i], aM_new_stars_r[i], 
					aM_baryons_r[i], aM_dark_r[i], aM_total_r[i], aM_total, 
					vmax, rmax, avg_gas_metallicity_II_r[i], 
					avg_gas_metallicity_Ia_r[i], avg_star_metallicity_II_r[i], 
					avg_star_metallicity_Ia_r[i], avg_new_star_metallicity_II_r[i], 
					avg_new_star_metallicity_Ia_r[i], avg_star_age_r[i] );

			fflush(rlist[i]);
		}

#ifdef HYDRO
		fprintf( bszlist, "%u %e %e %e %e", halos->list[ihalo].id, aM_gas, aM_dark, aM_total, Ysz );
		for ( i = 0; i < num_radii; i++ ) {
			fprintf( bszlist, " %e %e %e %e", aM_gas_r[i], aM_dark_r[i], aM_total_r[i], Ysz_r[i] );
		}
		fprintf( bszlist, "\n" );
		fflush(bszlist);
#endif /* HYDRO */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		fprintf( beszlist, "%u %e %e", halos->list[ihalo].id, Ysz, Yesz );
		for ( i = 0; i < num_radii; i++ ) {
			fprintf( beszlist, " %e %e", Ysz_r[i], Yesz_r[i] );
		}
		fprintf( beszlist, "\n" );
#endif

#ifdef ANALYSIS_XRAY
#ifdef ELECTRON_ION_NONEQUILIBRIUM 
		fprintf( btxlist, "%u %e %e", halos->list[ihalo].id, Tx, Tex );
		for ( i = 0; i < num_radii; i++ ) {
			fprintf( btxlist, " %e %e", Tx_r[i], Tex_r[i] );
		}
		fprintf( btxlist, "\n" );
#else
		fprintf( btxlist, "%u %e %e %e %e", halos->list[ihalo].id, aM_gas, aM_dark, aM_total, Tx );
		for ( i = 0; i < num_radii; i++ ) {
			fprintf( btxlist, " %e %e %e %e", aM_gas_r[i], aM_dark_r[i], aM_total_r[i], Tx_r[i] );
		}
		fprintf( btxlist, "\n" );
#endif
#endif /* ANALYSIS_XRAY */


		fprintf( bmpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id, rdout, rvdout, 
				aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir, vmax, rmax,
				avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia, 
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, 
				avg_star_age );

		/* write out profiles */
		for ( bin = 0; bin < num_bins; bin++ ) {
			fprintf( bmpro, "%.3f %.3f %e %e %e %e %e\n",
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					bin_total_dark_mass[bin]*aM0*hubble,
					bin_total_gas_mass[bin]*aM0*hubble,
					bin_total_cold_gas_mass[bin]*aM0*hubble,
					bin_total_star_mass[bin]*aM0*hubble,
					bin_total_new_star_mass[bin]*aM0*hubble );
		}

		fflush(bmpro);

		fprintf( bvpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id, 
				rdout, rvdout, aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir, vmax, rmax,
				avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia, 
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

		for ( bin = 0; bin < num_bins; bin++ ) {
			rcirc = rr[bin]*1000.0*r0/hubble; 
			vcirc_dm = 2.07498e-3 * sqrt( aM0*bin_total_dark_mass[bin] / rcirc );
			vcirc_gas = 2.07498e-3 * sqrt( aM0*bin_total_gas_mass[bin] / rcirc );
			vcirc_stars = 2.07498e-3 * sqrt( aM0*bin_total_star_mass[bin] / rcirc );
			vcirc_total = 2.07498e-3 * sqrt( aM0*(bin_total_dark_mass[bin]+
						bin_total_gas_mass[bin]+bin_total_star_mass[bin])/rcirc );

			profile->bin_dark_vrms[bin] = vfact*sqrt(profile->bin_dark_vrms[bin]);
			profile->bin_gas_vrms[bin] = vfact*sqrt(profile->bin_gas_vrms[bin]);
			profile->bin_star_vrms[bin] = vfact*sqrt(profile->bin_star_vrms[bin]);
			profile->bin_new_star_vrms[bin] = vfact*sqrt(profile->bin_new_star_vrms[bin]);

			fprintf( bvpro, "%.3f %.3f %e %e %e %e %.3f %.3f %.3f %.3f\n", 
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					profile->bin_dark_vrms[bin], 
					profile->bin_gas_vrms[bin], 
					profile->bin_star_vrms[bin], 
					profile->bin_new_star_vrms[bin],
					vcirc_dm, vcirc_gas, vcirc_stars, vcirc_total );
		}

		fflush(bvpro);

#ifdef HYDRO
		fprintf( bgpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id,
				rdout, rvdout, aM_gas, aM_cold_gas,
				aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir,
				vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia,
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

		for ( bin = 0; bin < num_bins; bin++ ) {
			fprintf( bgpro, "%.3f %.3f %e %e %e %e %e\n",
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					bin_total_gas_mass[bin]*aM0*hubble,
					bin_total_cold_gas_mass[bin]*aM0*hubble,
					profile->bin_gas_temperature[bin],
					profile->bin_gas_pressure[bin],
					profile->bin_gas_entropy[bin]
			       );
		}

		fflush(bgpro);

#ifdef MASS_WEIGHTED 
		fprintf( bgpro_massweighted, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id,
				rdout, rvdout, aM_gas, aM_cold_gas,
				aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir,
				vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia,
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

		for ( bin = 0; bin < num_bins; bin++ ) {
			fprintf( bgpro_massweighted, "%.3f %.3f %e %e %e %e %e\n",
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					bin_total_gas_mass[bin]*aM0*hubble,
					bin_total_cold_gas_mass[bin]*aM0*hubble,
					profile->bin_massweighted_temperature[bin],
					profile->bin_massweighted_pressure[bin],
					profile->bin_massweighted_entropy[bin]
				   );
		}

		fflush(bgpro_massweighted);
#endif /* MASS_WEIGHTED */


#ifdef COOLING
		fprintf( bcpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id,
				rdout, rvdout, aM_gas, aM_cold_gas,
				aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir,
				vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia,
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

		for ( bin = 0; bin < num_bins; bin++ ) {

    fprintf( bcpro, "# rmid rr [/h kpc] Mg Mcg [/h Msolar] dE_cool [ergs s^-1] <tcool> [Myr] <Zg_II> <Zg_Ia>\n");
			fprintf( bcpro, "%.3f %.3f %e %e %.3f %.3f %e %e\n",
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					bin_total_gas_mass[bin]*aM0*hubble,
					bin_total_cold_gas_mass[bin]*aM0*hubble,
					profile->bin_gas_coolingrate[bin],
					profile->bin_gas_tcool[bin],
					profile->bin_gas_metallicity_II[bin]/Zsolar,
					profile->bin_gas_metallicity_Ia[bin]/Zsolar
				   );
		}

		fflush(bgpro);
#endif /* COOLING */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		fprintf( bepro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id,
				rdout, rvdout, aM_gas, aM_cold_gas,
				aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir,
				vmax, rmax, avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia,
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, avg_star_age );

		for ( bin = 0; bin < num_bins; bin++ ) {
			fprintf( bepro, "%.3f %.3f %e %e %e %e %e %e %e %e %e\n",
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					bin_total_gas_mass[bin]*aM0*hubble,
					profile->bin_gas_temperature[bin],
					profile->bin_electron_temperature[bin],
					bin_total_sz_flux[bin],
					bin_total_electron_sz_flux[bin],
					profile->bin_electron_density[bin],
					profile->bin_electron_tei[bin],
					profile->bin_gas_pressure[bin],
					profile->bin_electron_pressure[bin]
			       );
		}

		fflush(bepro);
#endif /* ELECTRON_ION_NONEQULIBRIUM */
#endif /* HYDRO */

#ifdef STARFORM
		fprintf( bzpro, "# %u %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %.3f\n",
				halos->list[ihalo].id, rdout, rvdout, 
				aM_gas, aM_cold_gas, aM_stars, aM_new_stars, aM_baryons,
				aM_dark, aM_total, halos->list[ihalo].mvir, vmax, rmax,
				avg_gas_metallicity_II, avg_gas_metallicity_Ia,
				avg_star_metallicity_II, avg_star_metallicity_Ia, 
				avg_new_star_metallicity_II, avg_new_star_metallicity_Ia, 
				avg_star_age );

		for ( bin = 0; bin < num_bins; bin++ ) {
			if ( profile->bin_star_age[bin] > 0 ) {
				age_star = age_t( profile->bin_star_age[bin] );
			} else {
				age_star = -1.0;
			}

			fprintf( bzpro, "%.3f %.3f %e %e %e %e %e %e %.3f\n",
					rmid[bin]*rfact*hubble, rr[bin]*rfact*hubble,
					profile->bin_gas_metallicity_II[bin]/Zsolar, 
					profile->bin_gas_metallicity_Ia[bin]/Zsolar,
					profile->bin_star_metallicity_II[bin]/Zsolar, 
					profile->bin_star_metallicity_Ia[bin]/Zsolar,
					profile->bin_new_star_metallicity_II[bin]/Zsolar, 
					profile->bin_new_star_metallicity_Ia[bin]/Zsolar,
					age_star );
		}

		fflush(bzpro);
#endif /* STARFORM */

	}

        /* finish up */
        fclose( blist );
        fclose( bmpro );
        fclose( bvpro );

#ifdef HYDRO
        fclose( bszlist );
#ifdef ANALYSIS_XRAY
        fclose( btxlist );
#endif /* ANALYSIS_XRAY */
        fclose( bgpro );
#endif /* HYDRO */

#ifdef STARFORM
        fclose( bzpro );
#endif /* STARFORM */

#ifdef COOLING
		fclose( bcpro );
#endif /* COOLING */

#ifdef MASS_WEIGHTED 
		fclose( bgpro_massweighted );
#endif /* MASS_WEIGHTED */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		fclose( beszlist );
		fclose( bepro );
#endif

		for ( i = 0; i < num_radii; i++ ) {
			fclose( rlist[i] );
		}
}
