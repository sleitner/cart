#ifndef __RECIPE_POISSONRF12_H__
#define __RECIPE_POISSONRF12_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION
void star_formation_poissonRF12( int level, int icell, double dt_eff, double sfr );
#endif /* STAR_FORMATION */
#endif
