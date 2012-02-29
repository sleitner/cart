#ifndef __DENSITY_H__
#define __DENSITY_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)

void config_init_density();
void config_verify_density();

void initialize_density( int level );
void assign_density( int level );

#ifdef HYDRO
void assign_hydro_density( int level );
#endif

#ifdef PARTICLES
void assign_particle_density( int level );
void assign_particle_density_smoothed( int level );
#endif

#endif

#endif
