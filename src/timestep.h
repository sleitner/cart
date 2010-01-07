#ifndef __TIMESTEP_H__
#define __TIMESTEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>

extern int max_steps;
extern double timelimit;

void config_init_timestep();
void config_verify_timestep();

int global_timestep( double dt );
int timestep( int level, MPI_Comm local_comm );
void choose_timestep( double *dt );

#ifdef HYDRO
void hydro_timestep( int level, int *courant_cell, double *velocity );
#endif

extern int step;

extern double t_init;
extern double t_end;
DECLARE_LEVEL_ARRAY(double,dtl);
DECLARE_LEVEL_ARRAY(double,dtl_old);
DECLARE_LEVEL_ARRAY(double,tl);
DECLARE_LEVEL_ARRAY(double,tl_old);

#ifdef COSMOLOGY
extern double auni_init;
extern double auni_end;
DECLARE_LEVEL_ARRAY(double,abox);
DECLARE_LEVEL_ARRAY(double,abox_old);
DECLARE_LEVEL_ARRAY(double,auni);
#endif /* COSMOLOGY */

DECLARE_LEVEL_ARRAY(int,num_steps_on_level);

#endif
