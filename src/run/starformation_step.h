#ifndef __STAR_FORMATION_STEP_H__
#define __STAR_FORMATION_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

#ifdef HYDRO
void star_formation( int level, int time_multiplier );
#endif /* HYDRO */

void remap_star_ids();

#endif /* STAR_FORMATION */

#endif
