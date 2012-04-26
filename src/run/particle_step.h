#ifndef __PARTICLE_STEP_H__
#define __PARTICLE_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef PARTICLES 

#ifdef GRAVITY
void accelerate_particles( int level );
#endif
void move_particles( int level );

#endif /* PARTICLES */

#endif /* __PARTICLE_STEP_H__ */
