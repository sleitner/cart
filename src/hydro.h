#ifndef __HYDRO_H__
#define __HYDRO_H__

#include "defs.h"
#include "tree.h"

#ifdef HYDRO

#define gamma			(5.0/3.0)
#define gamma_max		10.0

#define COPY		0
#define RESTORE		1
#define COPY_ZERO_REF	2

#define	COPY_ALL_LEAFS		0
#define	COPY_SPLIT_NEIGHBORS	1
#define COPY_NO_SPLIT_NEIGHBORS	2

#ifdef PRESSURE_FLOOR
extern int pressure_floor_min_level;
extern float pressure_floor_factor;
#endif

extern float ref[num_cells];
extern int level_sweep_dir[max_level-min_level+1];

extern float gas_density_floor;
extern float gas_temperature_floor;

void hydro_step( int level );
void hydro_copy_vars( int level, int direction, int cell_type );
void apply_hydro_fluxes( int icell, double factor, double dxi_factor, double f[num_hydro_vars-1] );
void hydro_sweep_1d( int level );
#ifdef GRAVITY
void hydro_apply_gravity( int level );
#endif /* GRAVITY */
void compute_hydro_fluxes( int cell_list[4], double f[num_hydro_vars-1] );
void hydro_eos(int level);
void hydro_magic(int level);
void hydro_advance_internalenergy(int level);
void hydro_split_update( int level );

#endif

#endif
