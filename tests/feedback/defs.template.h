/*#define GCC_COMPILER*/
#define USER_PLUGIN

#define STARFORM 
#define PARTICLES 
#define COSMOLOGY  

#define HYDRO
#define REFINEMENT
/* #define MOMENTUM_DIFFUSION /\*affects refinement*\/ */
#define STAR_PRESSURE

#define num_root_grid_refinements	6
#define num_refinement_levels		3
#define num_octs			2000000	/* suitable for no refinement */
#define num_particles		        8	/* suitable for no refinement */
#define num_star_particles		8	/* suitable for no refinement */

#define ENRICH
#define ENRICH_SNIa

/* #define GRAVITY   */
#define COOLING
/* #define RADIATIVE_TRANSFER*/
