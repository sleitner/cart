#ifndef __STARFORMATION_H__
#define __STARFORMATION_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

#include "particle.h"

#ifdef STAR_PARTICLE_TYPES
extern int star_particle_type[num_star_particles];
#endif /* STAR_PARTICLE_TYPES */

#define STAR_TYPE_DELETED       (-1)
#define STAR_TYPE_NORMAL        0
#define STAR_TYPE_AGN           1
#define STAR_TYPE_STARII        2
#define STAR_TYPE_FAST_GROWTH   3


extern float star_formation_volume_min[nDim];
extern float star_formation_volume_max[nDim];

DECLARE_LEVEL_ARRAY(int,star_formation_frequency);

extern int num_local_star_particles;
extern particleid_t last_star_id;
extern int num_new_stars;

extern double total_stellar_mass;
extern double total_stellar_initial_mass;

extern float star_tbirth[num_star_particles];
extern float star_initial_mass[num_star_particles];

#ifdef ENRICHMENT
extern float star_metallicity_II[num_star_particles];
#ifdef ENRICHMENT_SNIa
extern float star_metallicity_Ia[num_star_particles];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

void config_init_star_formation();
void config_verify_star_formation();

void init_star_formation();
void star_formation_rate( int level, int num_level_cells, int *level_cells, float *sfr);
void grow_star_particle( int ipart, float dmass, int icell, int level);
int create_star_particle( int icell, float mass, double pdt, int type );

/* global parameters */

extern double sf_sampling_timescale;

#endif /* STAR_FORMATION */

#endif
