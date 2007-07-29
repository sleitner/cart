#ifndef __GRAVITY_H__
#define __GRAVITY_H__

#ifdef GRAVITY

void solve_poisson( int level, int flag);
void copy_potential( int level );
void interpolate_potential( int level ); 
void relax( int level, int flag);
void prolongate( int level );
void smooth( int level );
void restrict_to_level( int level );
void potential();

#ifdef HYDRO
void compute_accelerations_hydro( int level );
#endif /* HYDRO */

#ifdef PARTICLES
void compute_accelerations_particles( int level );
#endif /* PARTICLES */

#endif /* GRAVITY */

#endif
