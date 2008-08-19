#ifndef __TIMESTEP_H__
#define __TIMESTEP_H__

#include "defs.h"
#include "tree.h"

#define min_courant_velocity    1e-6

extern double a_init;
extern double a_end;
extern double t_init;
extern double t_end;

extern int max_steps;
extern int output_frequency;
extern int restart_frequency;
extern int particle_output_frequency;
extern int grid_output_frequency;
extern int tracer_output_frequency;

extern double cfl;
extern double particle_cfl;
extern double max_time_inc;
extern double min_time_dec;
extern double max_da;
extern double max_dt;
extern double max_frac_da;

int global_timestep( double dt );
int timestep( int level, MPI_Comm local_comm );
void choose_timestep( double *dt );

#ifdef HYDRO
void hydro_timestep( int level, int *courant_cell, double *velocity );
#endif

extern int step;
extern double dtl[max_level-min_level+1];
extern double dtl_old[max_level-min_level+1];
extern double tl[max_level-min_level+1];
extern double tl_old[max_level-min_level+1];
extern double aexp[max_level-min_level+1];
extern double aexp_old[max_level-min_level+1];

extern int num_steps_on_level[max_level-min_level+1];

#endif
