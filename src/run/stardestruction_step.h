#ifndef __STAR_DESTRUCTION_STEP_H__
#define __STAR_DESTRUCTION_STEP_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION
void star_destruction(int level);
#endif /* STAR_FORMATION */

#endif
