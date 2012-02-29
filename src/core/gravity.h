#ifndef __GRAVITY_H__
#define __GRAVITY_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef GRAVITY

void config_init_gravity();
void config_verify_gravity();

void solve_poisson( int level, int flag);
void relax( int level, int flag);
void prolongate( int level );
void smooth( int level );
void restrict_to_level( int level );

#endif /* GRAVITY */

#endif
