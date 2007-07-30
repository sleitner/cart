#ifndef __DEFS_H__
#define __DEFS_H__

#define HYDRO
#define LAPIDUS
#define DENSGRADSMOOTH
/* #define PARTICLES */
/* #define GRAVITY */
#define REFINEMENT
/* #define COSMOLOGY */
#define MOMENTUM_DIFFUSION
/* #define TEST_RIEMANN_SOLVER */
/* #define GRAVITY_IN_RIEMANN */
/* #define NDEBUG */
#define HYDRO_TRACERS

#define num_root_grid_refinements	6
#define num_refinement_levels		4
#define num_octs			500000	/* none needed with no refinement */
#define max_local_root_cells		(num_root_cells)
#define num_particles			0
#define SFC				HILBERT
#define T_min				(0.0)
#define num_tracers			(300000)

#endif
