#ifndef __FORMSTARS_POISSONRF12_H__
#define __FORMSTARS_POISSONRF12_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

void poissonRF12_config_init();
void poissonRF12_config_verify();
void poissonRF12_init();
void poissonRF12_setup(int level);
void poissonRF12_star_formation( int level, int icell, double dtl, double dt_eff, float sfr);

#endif /* STAR_FORMATION */
#endif
