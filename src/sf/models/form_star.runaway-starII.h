#ifndef __FORMSTARS_RUNAWAY_STARII_H__
#define __FORMSTARS_RUNAWAY_STARII_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION
void starII_runaway_config_init();
void starII_runaway_config_verify();
int starII_runaway_velocity(double mass_code, double vadd[nDim]);
#endif /* STAR_FORMATION */
#endif
