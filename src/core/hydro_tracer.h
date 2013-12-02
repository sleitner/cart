#ifndef __HYDRO_TRACER_H__
#define __HYDRO_TRACER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#if defined(HYDRO) && defined(HYDRO_TRACERS)

#include <stdint.h>
#include <limits.h>

#ifdef OLDSTYLE_32BIT_TRACERID
#define tracerid_t      int
#define MPI_TRACERID_T  MPI_INT
#define NULL_TRACER     (-1)
#define TRACERID_MAX    INT_MAX
#else
#define tracerid_t      int64_t
#define MPI_TRACERID_T  MPI_LONG
#define NULL_TRACER	    (-1L)
#define TRACERID_MAX    INT64_MAX
#endif

extern tracerid_t tracer_id[num_tracers];
extern double tracer_x[num_tracers][nDim];
extern int tracer_list_next[num_tracers];
extern int tracer_list_prev[num_tracers];

extern int cell_tracer_list[num_cells];

extern int num_tracer_row;
extern int num_local_tracers;
extern tracerid_t num_tracers_total;
extern int next_free_tracer;
extern int free_tracer_list;
extern int tracer_list_enabled;

extern int num_hydro_vars_traced;
extern int hydro_vars_traced[];
extern char *hydro_vars_traced_labels[];

int tracer_alloc( tracerid_t id );
void tracer_free( int tracer );
void tracer_list_free( int ihead );

void init_hydro_tracers();
void set_hydro_tracers( int min_tracer_level );

#ifdef PARTICLES
void set_hydro_tracers_to_particles();
#endif /* PARTICLES */

void update_tracer_list( int level );
void trade_tracer_lists( int *num_tracers_to_send, int *tracer_list_to_send, int trade_level );
void build_tracer_list();
void split_tracer_list( int icell );
void join_tracer_list( int icell );
void insert_tracer( int icell, int part );
void delete_tracer( int icell, int part );

#endif /* defined(HYDRO) && defined(HYDRO_TRACERS) */

#endif
