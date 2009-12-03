#ifndef __DENSITY_H__
#define __DENSITY_H__

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)

void initialize_density( int level );
void assign_density( int level );

#ifdef HYDRO
void assign_hydro_density( int level );
#endif

#ifdef PARTICLES
void assign_particle_density( int level );
#ifdef MAX_LEVEL_DARK_DENSITY
void assign_particle_density_smoothed( int level );
#endif
#endif

#endif

#endif
