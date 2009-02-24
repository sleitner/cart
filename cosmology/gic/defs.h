#ifndef __DEFS_H__
#define __DEFS_H__


#define DEBUG 2
#define DEBUG_MEMORY_USE


/* #define HYDRO */
#define GRAVITY
#define COOLING
/* #define STARFORM */
#define COSMOLOGY
#define PARTICLES
#define REFINEMENT
/* #define RADIATIVE_TRANSFER */
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

#define SF_RECIPE                       3


#define OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION 
#define OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID 


#define num_root_grid_refinements	5
#define num_refinement_levels		5
#define num_octs			100000
#define num_particles		        200000
#define num_star_particles              100000

#endif
