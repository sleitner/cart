#ifndef __DEFS_H__
#define __DEFS_H__


#define HYDRO
#define GRAVITY
#define COOLING
#define STARFORM
#define COSMOLOGY
#define PARTICLES
#define REFINEMENT
#define RADIATIVE_TRANSFER
/* #define CLOUDY_COOLING */

#define LAPIDUS
#define DENSGRADSMOOTH
#define PRESSURE_FLOOR
#define ADVECT_SPECIES

#define METALCOOLING
#define ENRICH
#define ENRICH_SNIa
#define FEEDBACK
#define FEEDBACK_SNIa
#define STELLARMASSLOSS
#define MinL_Jeans                      8
#define T_max_feedback                  2e7

#define SF_RECIPE                       3


#define OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION
#define OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID


#define num_root_grid_refinements	6
#define num_refinement_levels		9
#define num_octs			2000000
#define num_particles		        5000000
#define num_star_particles              1000000

#endif
