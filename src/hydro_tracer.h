#ifndef __HYDRO_TRACER_H__
#define __HYDRO_TRACER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(HYDRO_TRACERS)

#define NULL_TRACER (-1)

extern double tracer_x[num_tracers][nDim];
extern int tracer_id[num_tracers];
extern int tracer_list_next[num_tracers];
extern int tracer_list_prev[num_tracers];

extern int CELL_ARRAY(cell_tracer_list);

extern int num_tracer_row;
extern int num_local_tracers;
extern int num_tracers_total;
extern int next_free_tracer;
extern int free_tracer_list;
extern int tracer_list_enabled;

extern int num_hydro_vars_traced;
extern int hydro_vars_traced[];
extern char *hydro_vars_traced_labels[];

int tracer_alloc( int id );
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
