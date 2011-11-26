#include "config.h"

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "times.h"
#include "units.h"


double t_init = -1.0e38;
double t_end = 1.0e38;

DEFINE_LEVEL_ARRAY(double,tl);

/*
//  Time refinement scheme defaults to the old factor-of-two refinement.
*/
int min_time_refinement_factor = 2;
int max_time_refinement_factor = 2;
int time_refinement_level = min_level;

#ifdef COSMOLOGY
double auni_init = 1.0e-3;
double auni_end = 1.0;
DEFINE_LEVEL_ARRAY(double,abox);
DEFINE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

int max_steps = 0;
double timelimit = 0.0;

int max_mpi_sync_level = min_level;

#ifdef HYDRO 
double cfl_run = 0.6;
double cfl_max = 0.6;     /* max allowed, re-do the timestep if that number is exceeded. */
#endif /* HYDRO */

#ifdef PARTICLES
double particle_cfl = 0.0;
#endif /* PARTICLES */

double max_time_inc = 1.1;
double min_time_dec = 1.25;
double max_dt = 0.125;
#ifdef COSMOLOGY
double max_a_inc = 1.1;
double max_da = 3e-3;
#endif /* COSMOLOGY */

int step = 0;
int steps_before_increasing = 4;


/*
// NG: it is not clear how make all these 4 parameters consistent.
// For now keep them completely independent
//
#ifdef COSMOLOGY
void control_parameter_set_aini(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  t_init = tcode_from_auni(auni_init);
}


void control_parameter_set_aend(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  t_end = tcode_from_auni(auni_end);
}
#endif


void control_parameter_set_tini(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
#ifdef COSMOLOGY
  auni_init = auni_from_tcode(t_init);
#endif
}


void control_parameter_set_tend(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
#ifdef COSMOLOGY
  auni_end = auni_from_tcode(t_end);
#endif
}
*/


#ifdef HYDRO 
void control_parameter_set_cflrun(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  if(cfl_max < cfl_run) cfl_max = cfl_run;
}


void control_parameter_set_cflmax(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  if(cfl_max < cfl_run) cfl_max = cfl_run;
}
#endif /* HYDRO  */


#ifdef COSMOLOGY
void control_parameter_set_a_inc(const char *value, void *ptr, int ind)
{
  control_parameter_set_double(value,ptr,ind);
  if(ind == 2)
    {
      max_a_inc += 1.0; /* backward compatibility */
    }
}
#endif /* COSMOLOGY */


void config_init_timestep()
{
#ifdef COSMOLOGY
  ControlParameterOps control_parameter_aini = { control_parameter_set_double, control_parameter_list_double };
  ControlParameterOps control_parameter_aend = { control_parameter_set_double, control_parameter_list_double };
  ControlParameterOps control_parameter_a_inc = { control_parameter_set_a_inc, control_parameter_list_double };
#endif /* COSMOLOGY */
  ControlParameterOps control_parameter_tini = { control_parameter_set_double, control_parameter_list_double };
  ControlParameterOps control_parameter_tend = { control_parameter_set_double, control_parameter_list_double };
#ifdef HYDRO 
  ControlParameterOps control_parameter_cflrun = { control_parameter_set_cflrun, control_parameter_list_double };
  ControlParameterOps control_parameter_cflmax = { control_parameter_set_cflmax, control_parameter_list_double };
#endif /* HYDRO  */

#ifdef COSMOLOGY
  control_parameter_add3(control_parameter_aini,&auni_init,"auni-start","auni_init","a_init","starting value for the cosmic scale factor.");

  control_parameter_add3(control_parameter_aend,&auni_end,"auni-stop","auni_end","a_end","last value for the cosmic scale factor. The simulation stops if this value is reached.");

  control_parameter_add3(control_parameter_a_inc,&max_a_inc,"max-a-increment","max_a_inc","max_frac_da","the largest factor by which the cosmic scale factor is allowed to increase in one time-step.");

  control_parameter_add2(control_parameter_double,&max_da,"max-da","max_da","maximum allowed step in the cosmic scale factor.");
#endif /* COSMOLOGY */

  control_parameter_add2(control_parameter_tini,&t_init,"time-start","t_init","starting value for the code time variable (in code units).");

  control_parameter_add2(control_parameter_tend,&t_end,"time-stop","t_end","last value for the code time variable (in code units). The simulation stops if this value is reached.");

  control_parameter_add2(control_parameter_int,&max_steps,"num-steps","max_steps","number of time-steps to make. Zero value disables this limit.");

  control_parameter_add2(control_parameter_double,&timelimit,"walltime-limit","timelimit","the limit for the wall-clock time for the current job. Zero value disables this limit.");

  control_parameter_add3(control_parameter_int,&max_mpi_sync_level,"max-mpi-sync-level","max_mpi_sync_level","max_cfl_sync_level","maximum level at which MPI calls are synchronized across separate nodes.");

#ifdef HYDRO 
  control_parameter_add3(control_parameter_cflrun,&cfl_run,"cfl-run","cfl_run","cfl","the CFL number for setting the time-step.");

  control_parameter_add3(control_parameter_cflmax,&cfl_max,"cfl-max","cfl_max","cfl","the maximum acceptable CFL number. If this number is exceeded, the time-step needs to be redone. It is a good sense to set this number just a little bit higher than <cfl-run>, to avoid extra restarts due to numerical noise.");
#endif /* HYDRO */

#ifdef PARTICLES
  control_parameter_add2(control_parameter_double,&particle_cfl,"particle-cfl","particle_cfl","the CFL number for particle dynamics. In HYDRO mode this number is usually not needed, as the grid CFL conditions superceeds that of particles. Setting it to zero disables this limit.");
#endif /* PARTICLES */

  control_parameter_add2(control_parameter_double,&max_time_inc,"max-timestep-increment","max_time_inc","the largest factor by which the time-step is allowed to increase.");

  control_parameter_add2(control_parameter_double,&min_time_dec,"min-timestep-decrement","min_time_dec","the smallest factor by which the time-step is allowed to decrease.");

  control_parameter_add2(control_parameter_double,&max_dt,"max-dt","max_dt","maximum allowed time-step.");

  control_parameter_add2(control_parameter_int,&steps_before_increasing,"timesteps-before-increasing","steps_before_increasing","number of global time-steps to make before the time-step is allowed to increase.");

  control_parameter_add2(control_parameter_int,&min_time_refinement_factor,"time-refinement-factor:min","min_time_refinement_factor","minimum allowed refinement factor between the time-steps on two successive levels.");

  control_parameter_add2(control_parameter_int,&max_time_refinement_factor,"time-refinement-factor:max","max_time_refinement_factor","maximum allowed refinement factor between the time-steps on two successive levels. If <time-refinement-factor:min> = <time-refinement-factor:max>, then time refinement is done by a constant factor; <time-refinement-factor:min> = <time-refinement-factor:max> = 2 recovers the old factor-of-two scheme.");

  control_parameter_add2(control_parameter_int,&time_refinement_level,"time-refinement-level","time_refinement_level","lowest level at which time refinement is allowed; all higher levels take the same time-steps as this level (a-la FLASH).");
}


void config_verify_timestep()
{
#ifdef COSMOLOGY
  cart_assert(auni_init>0.0 && !(auni_init>auni_end));

  cart_assert(max_a_inc > 1.0);

  cart_assert(max_da >= 0.0);
#endif /* COSMOLOGY */

  cart_assert(!(t_init > t_end));

  cart_assert(max_steps >= 0);

  cart_assert(max_mpi_sync_level >= min_level);

#ifdef HYDRO 
  cart_assert(cfl_run>0.0 && !(cfl_run>cfl_max));
#endif /* HYDRO */

#ifdef PARTICLES
  cart_assert(!(particle_cfl < 0.0));
#endif /* PARTICLES */

  cart_assert(max_time_inc > 1.0);

  cart_assert(min_time_dec > 1.0);

  cart_assert(max_dt >= 0.0);

  cart_assert(steps_before_increasing > 0);

  cart_assert(min_time_refinement_factor > 0);
  cart_assert(max_time_refinement_factor >= min_time_refinement_factor);
  cart_assert(time_refinement_level >= min_level);

#ifndef HYDRO
  if(min_time_refinement_factor==1 && max_time_refinement_factor>1 && particle_cfl==0.0)
    {
      cart_error("In a pure N-body mode, with <particle-cfl> parameter set to 0.0, it is incorrect to set the <time-refinement-factor:min> parameter to 1. In that case there is no condition to make time-steps on lower levels smaller than on the top level, so the particle trajectories will be integrated incorrectly.");
    }
#endif /* HYDRO */
}

