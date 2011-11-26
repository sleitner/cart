#ifndef __STARFORMATION_STEP_H__
#define __STARFORMATION_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

#ifdef HYDRO
void star_formation( int level, int time_multiplier );
void create_star_particle( int icell, float mass );
#endif /* HYDRO */

void remap_star_ids();

#endif /* STARFORM */

#endif
