#ifndef __DEFS_H__
#define __DEFS_H__

#define PARTICLES
#define GRAVITY
#define COSMOLOGY
#define HYDRO

#define REFINEMENT

/* Hydrodynamic parameters, should be on by default */
#define LAPIDUS
#define DENSGRADSMOOTH
#define PRESSURE_FLOOR

#define num_root_grid_refinements	7
#define num_refinement_levels		6
#define num_octs			100000
#define num_particles			(2000000)

#endif