#ifndef __HYDRO_H__
#define __HYDRO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef HYDRO

#define COPY		0
#define RESTORE		1
#define COPY_ZERO_REF	2

#define	COPY_ALL_LEAFS		0
#define	COPY_SPLIT_NEIGHBORS	1
#define COPY_NO_SPLIT_NEIGHBORS	2

extern int pressure_floor_min_level;

extern float ref[num_cells];
DECLARE_LEVEL_ARRAY(int,level_sweep_dir);

extern float gas_density_floor;
extern float gas_temperature_floor;

void config_init_hydro();
void config_verify_hydro();

void hydro_step( int level );
void hydro_copy_vars( int level, int direction, int cell_type );
void apply_hydro_fluxes( int icell, double factor, double dxi_factor, double f[ /* num_hydro_vars-1 */ ] );
void hydro_sweep_1d( int level );
#ifdef GRAVITY
void hydro_apply_gravity( int level );
#endif /* GRAVITY */
void compute_hydro_fluxes( int cell_list[4], double f[ /* num_hydro_vars-1 */ ] );
void hydro_eos(int level);
void hydro_magic_one_cell(int level);
void hydro_magic(int level);
void hydro_advance_internalenergy(int level);
void hydro_split_update( int level );

float cell_gas_kinetic_energy(int cell);
float cell_gas_temperature(int cell);

#endif /* HYDRO */
#endif /* __HYDRO_H__ */
