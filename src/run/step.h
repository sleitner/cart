#ifndef __STEP_H__
#define __STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


void run(int restart, const char *restart_label);
void set_timestepping_scheme();


DECLARE_LEVEL_ARRAY(double,dtl);
DECLARE_LEVEL_ARRAY(double,dtl_old);

DECLARE_LEVEL_ARRAY(int,time_refinement_factor);
DECLARE_LEVEL_ARRAY(int,time_refinement_factor_old);

DECLARE_LEVEL_ARRAY(unsigned int,num_steps_on_level);

DECLARE_LEVEL_ARRAY(double,tl_old);
#ifdef COSMOLOGY
DECLARE_LEVEL_ARRAY(double,abox_old);
#endif /* COSMOLOGY */

extern int step;

#endif
