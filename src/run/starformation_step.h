#ifndef __STARFORMATION_STEP_H__
#define __STARFORMATION_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

#ifdef HYDRO
void star_formation( int level, int time_multiplier );
#endif /* HYDRO */

void remap_star_ids();

#endif /* STARFORM */

#endif
