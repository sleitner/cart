#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "defs.h"
#include "auxiliary.h"
#include "halos.h"

#define subhalo_exclude_radius	(1.5)
#define EXCLUDE_SUBHALOS

#define MASS_WEIGHTED 

/*
#define rbinmax         (4e3/hubble)
#define rbinmin         (40./hubble)
#define max_bins        20
*/

#define rbinmax         (2e4/hubble)
#define rbinmin         (2.0/hubble)
#define max_bins        100

#define points_per_cell 100

#define Tcold           (5e4)
#define new_star_age    (0.1)

typedef struct {
	double rlmin;
	double drl;
	int num_bins; 

	int bin_dark_num[max_bins];
	double bin_dark_mass[max_bins];
	double bin_dark_momentum[nDim][max_bins];
	double bin_dark_vrms[max_bins];

	int bin_star_num[max_bins];
	double bin_star_mass[max_bins];
	double bin_star_momentum[nDim][max_bins];
	double bin_star_vrms[max_bins];
	double bin_star_age[max_bins];
	double bin_star_metallicity_II[max_bins];
	double bin_star_metallicity_Ia[max_bins];

	int bin_new_star_num[max_bins];
	double bin_new_star_mass[max_bins];
	double bin_new_star_momentum[nDim][max_bins];
	double bin_new_star_vrms[max_bins];
	double bin_new_star_metallicity_II[max_bins];
	double bin_new_star_metallicity_Ia[max_bins];

#ifdef HYDRO
	double bin_gas_mass[max_bins];
	double bin_volume[max_bins];
	double bin_excluded_gas_mass[max_bins];
	double bin_excluded_volume[max_bins];
	double bin_gas_velocity[nDim][max_bins];
	double bin_gas_vrms[max_bins];
	double bin_cold_gas_mass[max_bins];
	double bin_gas_temperature[max_bins];
	double bin_gas_entropy[max_bins];
	double bin_gas_pressure[max_bins];
	double bin_gas_coolingrate[max_bins];
	double bin_gas_tcool[max_bins];
	double bin_gas_metallicity_II[max_bins];
	double bin_gas_metallicity_Ia[max_bins];
	double bin_sz_flux[max_bins];

#ifdef MASS_WEIGHTED
	double bin_massweighted_temperature[max_bins];
	double bin_massweighted_entropy[max_bins];
	double bin_massweighted_pressure[max_bins];
#endif /* MASS_WEIGHTED */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double bin_electron_pressure[max_bins];
	double bin_electron_temperature[max_bins];
	double bin_electron_sz_flux[max_bins];
	double bin_electron_tei[max_bins];
	double bin_electron_density[max_bins];
#endif
	
#ifdef ANALYSIS_XRAY
	double bin_xray_Fcont[max_bins];
	double bin_xray_Fline[max_bins];
	double bin_xray_avgE[max_bins];
	double bin_xray_Tcont1[max_bins];
	double bin_xray_Tcont2[max_bins];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double bin_exray_Fcont[max_bins];
	double bin_exray_Fline[max_bins];
	double bin_exray_avgE[max_bins];
	double bin_exray_Tcont1[max_bins];
	double bin_exray_Tcont2[max_bins];
#endif

	double bin_xray_wsl[max_bins];
	double bin_xray_wslt[max_bins];
	double bin_xray_Lxray[max_bins];
#endif /* ANALYSIS_XRAY */
#endif
} profile_struct;

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos, halo_list *subhalos );

#endif
