#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "defs.h"
#include "auxiliary.h"
#include "halos.h"

#define points_per_cell	10

#define rbinmin		0.05
#define rbinmax		1.0
#define max_bins	200 

#define phase_rhobins 		300
#define phase_Tbins		300
#define phase_Sbins		300

#define	phase_Tlmin	(-1.0)
#define	phase_Tlmax 	11.0
#define	phase_dlt	((phase_Tlmax-phase_Tlmin)/phase_Tbins)

#define	phase_rhoglmin (-10.0)
#define	phase_rhoglmax (2.0)
#define	phase_drhol	((phase_rhoglmax-phase_rhoglmin)/phase_rhobins)

#define phase_Slmin	-1.0
#define phase_Slmax	5.0
#define phase_dls	((phase_Slmax-phase_Slmin)/phase_Sbins)

typedef struct {
        double rlmin;
        double drl;
        int num_bins;

	float phase_table[phase_Tbins][phase_rhobins];
	float entropy_table[phase_Sbins][phase_rhobins];

	float radial_temperature[phase_Tbins][max_bins];
	float radial_entropy[phase_Sbins][max_bins];
	float radial_density[phase_rhobins][max_bins];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	float electron_phase_table[phase_Tbins][phase_rhobins];
	float radial_electron_temperature[phase_Tbins][max_bins];
#endif
} analysis_struct;

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos );

#endif
