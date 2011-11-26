#ifndef __HYDRO_H__
#define __HYDRO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef HYDRO

extern int pressure_floor_min_level;

extern float gas_density_floor;
extern float gas_temperature_floor;

void config_init_hydro();
void config_verify_hydro();

float cell_gas_kinetic_energy(int cell);
float cell_gas_temperature(int cell);
float cell_radiation_pressure(int cell);

float cell_sobolev_length2(int cell, int level, float *vel);
#define cell_sobolev_length(c) cell_sobolev_length2(c,cell_level(c),NULL)

#endif /* HYDRO */
#endif /* __HYDRO_H__ */
