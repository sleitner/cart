#ifndef __TIMESTEP_H__
#define __TIMESTEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


extern int max_steps;
extern double timelimit;

void config_init_timestep();
void config_verify_timestep();

int global_timestep();
void set_timestepping_scheme();

#ifdef HYDRO
void hydro_cfl_condition( int level, int *courant_cell, double *velocity );
#endif

extern int step;

extern double t_init;
extern double t_end;
DECLARE_LEVEL_ARRAY(double,dtl);
DECLARE_LEVEL_ARRAY(double,dtl_old);
DECLARE_LEVEL_ARRAY(double,tl);
DECLARE_LEVEL_ARRAY(double,tl_old);

DECLARE_LEVEL_ARRAY(int,time_refinement_factor);
DECLARE_LEVEL_ARRAY(int,time_refinement_factor_old);

#ifdef COSMOLOGY
extern double auni_init;
extern double auni_end;
DECLARE_LEVEL_ARRAY(double,abox);
DECLARE_LEVEL_ARRAY(double,abox_old);
DECLARE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

DECLARE_LEVEL_ARRAY(unsigned int,num_steps_on_level);

#endif
