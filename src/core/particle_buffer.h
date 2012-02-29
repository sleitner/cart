#ifndef __PARTICLE_BUFFER_H__
#define __PARTICLE_BUFFER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#include "starformation.h"

#ifdef PARTICLES 
#ifdef STAR_PARTICLE_TYPES
void build_particle_buffer( int specie, int star_particle_type );
#else
void build_particle_buffer( int specie );
#endif /* STAR_PARTICLE_TYPES */
void destroy_particle_buffer();
#endif /* PARTICLES */
#endif /* __PARTICLE_BUFFER_H__ */
