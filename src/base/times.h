#ifndef __TIMES_H__
#define __TIMES_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


DECLARE_LEVEL_ARRAY(double,tl);

#ifdef COSMOLOGY
DECLARE_LEVEL_ARRAY(double,abox);
DECLARE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

void config_init_times();
void config_verify_times();

#endif /* __TIMES_H__ */
