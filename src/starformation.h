#ifndef __STARFORMATION_H__
#define __STARFORMATION_H__

#include "defs.h"
#include "tree.h"


#ifdef STARFORM

extern float star_formation_volume_min[nDim];
extern float star_formation_volume_max[nDim];

extern int star_formation_frequency[max_level-min_level+1];

extern int num_local_star_particles;
extern int last_star_id;
extern int num_new_stars;

extern double total_stellar_mass;
extern double total_stellar_initial_mass;

extern float star_tbirth[num_star_particles];
extern float star_initial_mass[num_star_particles];

#ifdef ENRICH
extern float star_metallicity_II[num_star_particles];
#ifdef ENRICH_SNIa
extern float star_metallicity_Ia[num_star_particles];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

void init_star_formation();
void star_formation( int level, int time_multiplier );
void create_star_particle( int icell, float mass );
void remap_star_ids();

/* global parameters */
extern double C_fb;
extern double fmass_met;
extern double C_fbIa;
extern double RIaf;

extern double alpha_SF;
extern double eps_SF;
extern double dtmin_SF;
extern double tau_SF;
extern double dm_star_min;
extern double rho_SF;
extern double T_SF;
extern double a_IMF;
extern double aM_stl;
extern double aM_stu;
extern double aM_SNII;
extern double aM_SNIa1;
extern double aM_SNIa2;
extern double ejM_SNIa;
extern double C_SNIa;
extern double t_SNIa;
extern double E_51;
extern double t_fb;
extern double T0_ml;
extern double c0_ml;

extern double T_max_feedback;

#endif /* STARFORM */

#endif
