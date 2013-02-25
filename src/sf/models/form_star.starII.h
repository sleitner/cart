#ifndef __FORMSTARS_STARII_H__
#define __FORMSTARS_STARII_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION
void starII_config_init();
void starII_config_verify();
void starII_init();
void starII_setup(int level);
void starII_star_formation( int level, int icell, double dtl, double dt_eff, float sfr);

#endif /* STAR_FORMATION */
#endif
