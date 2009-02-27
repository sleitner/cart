#ifndef __DEFS_H__
#define __DEFS_H__

/* #define HYDRO  */
#define PARTICLES
#define GRAVITY
#define COSMOLOGY

/* #define PRESSURELESS_FLUID */
/* #define MOMENTUM_DIFFUSION */
/* #define TEST_RIEMANN_SOLVER */
/* #define GRAVITY_IN_RIEMANN */
/* #define NDEBUG */

#define num_root_grid_refinements	9
#define num_refinement_levels		0
#define num_octs			0
#define max_local_root_cells		(num_root_cells/6)
#define num_particles			(num_grid*num_grid*num_grid/6)

#endif