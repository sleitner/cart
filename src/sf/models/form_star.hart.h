#ifndef __FORMSTARS_HART_H__
#define __FORMSTARS_HART_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION
void formstar_hart_config_init();
void formstar_hart_config_verify();
void formstar_hart_setup(int level);

#endif /* STAR_FORMATION */
#endif
