#define RADIATIVE_TRANSFER 

#define STAR_FORMATION 
#define SF_RECIPE <hart>
#define PARTICLES 
#define COSMOLOGY  

#define HYDRO
#define REFINEMENT
#define PREFIX_JOBNAME_TO_OUTPUT_FILES


#define num_root_grid_refinements	4
#define num_refinement_levels		3
#define num_octs			2000000	/* suitable for no refinement */
#define num_particles		        100	/* suitable for no refinement */
#define num_star_particles		100	/* suitable for no refinement */

#define ENRICH
#define ENRICH_SNIa

/* #define GRAVITY   */
#define COOLING

/* #define STAR_PRESSURE */
/* #define STAR_PRESSURE_FROM_PARTICLES  */
/* /\* #define STAR_PRESSURE_FROM_PARTICLES_TO_CELLS *\/ */
/* #define STAR_PRESSURE_TO_VELOCITY  */

/* /\* #define STAR_PRESSURE_IN_INTERNAL_ENERGY *\/ */



/* #ifdef STAR_PRESSURE_FROM_PARTICLES */
/* #ifndef STAR_PRESSURE_FROM_PARTICLES_TO_CELLS */
/* #ifndef STAR_PRESSURE_TO_VELOCITY */
/* #error "define STAR_PRESSURE_FROM_PARTICLES_TO_CELLS in hydro" */
/* #endif */
/* #endif */
/* #endif */

/* #if defined(STAR_PRESSURE_FROM_PARTICLES) && defined(STAR_PRESSURE_IN_INTERNAL_ENERGY) */
/* #error "STAR_PRESSURE IE or Particles" */
/* #endif */
