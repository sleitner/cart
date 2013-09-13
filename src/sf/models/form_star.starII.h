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
void starII_creation( double dmstarII, int icell, int level, double dtl );

#endif /* STAR_FORMATION */
#endif
