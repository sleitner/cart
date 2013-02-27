#ifndef __CELL_BUFFER_H__
#define __CELL_BUFFER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "index_hash.h"
#include "oct_hash.h"

extern int root_buffer_enabled;
extern int buffer_enabled;
extern int *buffer_cell_hash_index;
extern int *buffer_cell_sfc_index;

extern index_hash *buffer_root_hash;
extern oct_hash *buffer_oct_hash[MAX_PROCS];
extern oct_hash *buffer_oct_reverse_hash[MAX_PROCS];

DECLARE_LEVEL_ARRAY(int,num_buffer_cells);
DECLARE_LEVEL_ARRAY(int,buffer_oct_list);

DECLARE_LEVEL_ARRAY(int*,num_remote_buffers);
DECLARE_LEVEL_ARRAY(int**,remote_buffers);
DECLARE_LEVEL_ARRAY(int*,num_local_buffers);
DECLARE_LEVEL_ARRAY(int**,local_buffers);

void init_cell_buffer();
int cell_buffer_hash( int sfc );
int cell_buffer_exists( int index );
int cell_buffer_local_index( int index );
void build_root_cell_buffer();
void build_cell_buffer();
void destroy_cell_buffer();

void update_buffer_level( int level, const int *var_indices, int num_update_vars );

#ifdef HYDRO
void merge_buffer_cell_gas_density_momentum( int level );
#endif
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
void merge_buffer_cell_densities( int level );
#endif

void split_buffer_cells( int level, int *cells_to_split, int num_cells_to_split );
void join_buffer_cells( int level, int *octs_to_join, int *parent_root_sfc, int num_octs_to_join );

#define PROC_ANY	-1

int cell_can_prune( int cell, int proc );

#endif
