#ifndef __FORMSTARS_CONTINUOUS_H__
#define __FORMSTARS_CONTINUOUS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION
void continuous_config_init();
void continuous_config_verify();
void continuous_init();
void continuous_setup(int level);
void continuous_star_formation( int level, int icell, double dtl, double dt_eff, float sfr);

#endif /* STAR_FORMATION */
#endif
