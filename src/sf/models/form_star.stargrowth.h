#ifndef __FORMSTARS_STARGROWTH_H__
#define __FORMSTARS_STARGROWTH_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef STAR_FORMATION
void grow_star_particle( int ipart, float delta_mass, int icell, int level);
#endif /* STAR_FORMATION */
#endif
