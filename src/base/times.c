#include "config.h"

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "times.h"
#include "units.h"


DEFINE_LEVEL_ARRAY(double,tl);
DEFINE_LEVEL_ARRAY(double,tl_old);
DEFINE_LEVEL_ARRAY(double,dtl);
DEFINE_LEVEL_ARRAY(double,dtl_old);

DEFINE_LEVEL_ARRAY(int,time_refinement_factor);
DEFINE_LEVEL_ARRAY(int,time_refinement_factor_old);

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
DEFINE_LEVEL_ARRAY(double,abox_old);
#else
double t_init = -1.0e38;
double t_end = 1.0e38;
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

double max_dt_inc = 1.1;
double min_dt_dec = 1.25;
double tol_dt_grow = 2.0;
double max_dt = 0.125;
#ifdef COSMOLOGY
double max_a_inc = 1.1;
double max_da = 3e-3;
#endif /* COSMOLOGY */

double reduce_dt_factor_shallow_dec = 0.2;
double reduce_dt_factor_deep_dec = 0.8;

int step = 0;
int current_step_level = -1;

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


void config_init_times()
{
  int level;

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
#else
  control_parameter_add2(control_parameter_tini,&t_init,"time-start","t_init","starting value for the code time variable (in code units).");

  control_parameter_add2(control_parameter_tend,&t_end,"time-stop","t_end","last value for the code time variable (in code units). The simulation stops if this value is reached.");
#endif /* COSMOLOGY */

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

  control_parameter_add3(control_parameter_double,&max_dt_inc,"max-timestep-increment","max_dt_inc","max_time_inc","the largest factor by which the time-step is allowed to increase.");

  control_parameter_add3(control_parameter_double,&min_dt_dec,"min-timestep-decrement","min_dt_dec","min_time_dec","the smallest factor by which the time-step is allowed to decrease.");

  control_parameter_add2(control_parameter_double,&tol_dt_grow,"tolerance-for-timestep-increase","tol_dt_grow","the tolerance factor by which the optimal time-step should increase the time-step actually taken on the previous step to be allowed to increase it value; if the time-step on a given level <dtl> is between <dtl_old> and <tolerance-for-timestep-increase>*<dtl_old>, then <dtl> is set to <dtl_old>.");

  control_parameter_add2(control_parameter_double,&max_dt,"max-dt","max_dt","maximum allowed time-step.");

  control_parameter_add2(control_parameter_double,&reduce_dt_factor_shallow_dec,"reduce-timestep-factor:shallow-decrement","reduce_dt_factor_shallow_dec","a factor by which, after the CFL violation, the time-step reduction factor is decreased for levels shallower than the offending level. For example, if the CFL violation requires a time-step reduction of a factor of 3 at level 5, timesteps on levels from 0 to 4 will be reduced by factors 1+2*<reduce-timestep-factor:shallow-decrement>^(5-level). Must be between 0 and 1 inclusive.");

  control_parameter_add2(control_parameter_double,&reduce_dt_factor_deep_dec,"reduce-timestep-factor:deep-decrement","reduce_dt_factor_deep_dec","a factor by which, after the CFL violation, the time-step reduction factor is decreased for levels deeper than the offending level. For example, if the CFL violation requires a time-step reduction of a factor of 3 at level 5, timesteps on levels 6 and deeper will be reduced by factors 1+2*<reduce-timestep-factor:deep-decrement>^(level-5). Must be between 0 and 1 inclusive.");

  control_parameter_add2(control_parameter_int,&min_time_refinement_factor,"time-refinement-factor:min","min_time_refinement_factor","minimum allowed refinement factor between the time-steps on two successive levels.");

  control_parameter_add2(control_parameter_int,&max_time_refinement_factor,"time-refinement-factor:max","max_time_refinement_factor","maximum allowed refinement factor between the time-steps on two successive levels. If <time-refinement-factor:min> = <time-refinement-factor:max>, then time refinement is done by a constant factor; <time-refinement-factor:min> = <time-refinement-factor:max> = 2 recovers the old factor-of-two scheme.");

  control_parameter_add2(control_parameter_int,&time_refinement_level,"time-refinement-level","time_refinement_level","lowest level at which time refinement is allowed; all higher levels take the same time-steps as this level (a-la FLASH).");

  for(level=min_level; level<=max_level; level++)
    {
      tl_old[level] = 0.0;
      dtl[level] = dtl_old[level] = 0.0;
      time_refinement_factor[level] = time_refinement_factor_old[level] = 1;
#ifdef COSMOLOGY
      abox_old[level] = 0.0;
#endif
    }
}


void config_verify_times()
{
#ifdef COSMOLOGY
  VERIFY(auni-stop, 1 );

  VERIFY(auni-start, auni_init>0.0 && !(auni_init>auni_end) );

  VERIFY(max-a-increment, max_a_inc > 1.0 );

  VERIFY(max-da, max_da >= 0.0 );
#else
  VERIFY(time-stop, 1 );
  VERIFY(time-start, !(t_init > t_end) );
#endif /* COSMOLOGY */

  VERIFY(walltime-limit, !(timelimit < 0.0) );

  VERIFY(num-steps, max_steps >= 0 );

  VERIFY(max-mpi-sync-level, max_mpi_sync_level >= min_level );

#ifdef HYDRO 
  VERIFY(cfl-max, cfl_max > 0.0 );
  VERIFY(cfl-run, cfl_run>0.0 && !(cfl_run>cfl_max) );
#endif /* HYDRO */

#ifdef PARTICLES
  VERIFY(particle-cfl, !(particle_cfl < 0.0) );
#endif /* PARTICLES */

  VERIFY(max-timestep-increment, max_dt_inc > 1.0 );

  VERIFY(min-timestep-decrement, min_dt_dec > 1.0 );

  VERIFY(tolerance-for-timestep-increase, !(tol_dt_grow < 1.0) );

  VERIFY(max-dt, max_dt > 0.0 ); /* used in the first cfl*/

  VERIFY(reduce-timestep-factor:shallow-decrement, reduce_dt_factor_shallow_dec>=0.0 && reduce_dt_factor_shallow_dec<=1.0 );
  VERIFY(reduce-timestep-factor:deep-decrement, reduce_dt_factor_deep_dec>=0.0 && reduce_dt_factor_deep_dec<=1.0 );

  VERIFY(time-refinement-factor:min, min_time_refinement_factor > 0 );
  VERIFY(time-refinement-factor:max, max_time_refinement_factor >= min_time_refinement_factor );
  VERIFY(time-refinement-level, time_refinement_level >= min_level );

#if !defined(HYDRO) && defined(PARTICLES)
  if(min_time_refinement_factor==1 && max_time_refinement_factor>1 && particle_cfl==0.0)
    {
      cart_error("In a pure N-body mode, with <particle-cfl> parameter set to 0.0, it is incorrect to set the <time-refinement-factor:min> parameter to 1. In that case there is no condition to make time-steps on lower levels smaller than on the top level, so the particle trajectories will be integrated incorrectly.");
    }
#endif /* !HYDRO && PARTICLES */
}

