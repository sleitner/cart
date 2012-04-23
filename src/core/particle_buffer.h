#ifndef __PARTICLE_BUFFER_H__
#define __PARTICLE_BUFFER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#include "starformation.h"

#ifdef PARTICLES 
void build_particle_buffer( int specie, int subspecies );
void destroy_particle_buffer();
#endif /* PARTICLES */
#endif /* __PARTICLE_BUFFER_H__ */
