#define HYDRO 
/* #define PARTICLES  */
#define COSMOLOGY

#define GRAVITY 
#define GRAVITY_IN_RIEMANN 

#ifdef DEBUG17
#define LAPIDUS 
#define DENSGRADSMOOTH 
#define ADVECT_SPECIES
#endif

/* #define MOMENTUM_DIFFUSION */
/* #define TEST_RIEMANN_SOLVER */
/* #define GRAVITY_IN_RIEMANN */
/* #define NDEBUG */

#define num_root_grid_refinements	6
#define num_refinement_levels		0
#define num_octs			2000000	
#define num_particles			(num_grid*num_grid*num_grid) 
