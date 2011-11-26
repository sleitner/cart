#ifndef __TIMES_H__
#define __TIMES_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


void config_init_timestep();
void config_verify_timestep();


extern int step;

DECLARE_LEVEL_ARRAY(double,tl);

#ifdef COSMOLOGY
extern double auni_init;
DECLARE_LEVEL_ARRAY(double,abox);
DECLARE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

extern double t_init;
extern double max_dt;

#endif /* __TIMES_H__ */
