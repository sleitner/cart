#ifndef __DEFS_H__
#define __DEFS_H__

#define HYDRO
#define LAPIDUS
#define DENSGRADSMOOTH

/*
#define HYDRO_TRACERS
#define num_tracers		(num_root_cells)
*/

#define num_root_grid_refinements	6
#define num_refinement_levels		0
#define num_octs			(num_root_cells/num_children)
#define SFC				SLAB

#endif
