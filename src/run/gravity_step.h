#ifndef __GRAVITY_STEP_H__
#define __GRAVITY_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef GRAVITY

#ifdef HYDRO
void copy_potential( int level );
void interpolate_potential( int level ); 
void compute_accelerations_hydro( int level );
#endif /* HYDRO */

#ifdef PARTICLES
void compute_accelerations_particles( int level );
#endif /* PARTICLES */

#endif /* GRAVITY */

#endif
