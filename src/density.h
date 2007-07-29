#ifndef __DENSITY_H__
#define __DENSITY_H__

#ifdef GRAVITY

void initialize_density( int level );
void assign_density( int level );

#ifdef HYDRO
void assign_hydro_density( int level );
#endif

#ifdef PARTICLES
void assign_particle_density( int level );
#endif

#endif

#endif
