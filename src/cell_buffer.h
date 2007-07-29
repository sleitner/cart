#ifndef __CELL_BUFFER_H__
#define __CELL_BUFFER_H__

#include "tree.h"
#include "index_hash.h"

extern int root_buffer_enabled;
extern int buffer_enabled;
extern int *buffer_cell_hash_index;
extern int *buffer_cell_sfc_index;

extern index_hash *buffer_root_hash;
extern index_hash *buffer_oct_hash[MAX_PROCS];
extern index_hash *buffer_oct_reverse_hash[MAX_PROCS];

extern int num_buffer_cells[max_level-min_level+1];
extern int buffer_oct_list[max_level-min_level+1]; 
extern int buffer_oct_list_needs_ordering[max_level-min_level+1];

extern int *remote_buffers[max_level-min_level+1][MAX_PROCS];
extern int num_remote_buffers[max_level-min_level+1][MAX_PROCS];
extern int num_local_buffers[max_level-min_level+1][MAX_PROCS];
extern int *local_buffers[max_level-min_level+1][MAX_PROCS];

void init_cell_buffer();
int cell_buffer_hash( int sfc );
int cell_buffer_exists( int index );
int cell_buffer_local_index( int index );
void build_root_cell_buffer();
void build_cell_buffer();
void destroy_cell_buffer();

void update_buffer_level( int level, const int *var_indices, int num_update_vars );

#ifdef GRAVITY
void merge_buffer_cell_densities( int level );
#endif

void split_buffer_cells( int level, int *cells_to_split, int num_cells_to_split );
void join_buffer_cells( int level, int *octs_to_join, int *parent_root_sfc, int num_octs_to_join );

#define PROC_ANY	-1

int cell_can_prune( int cell, int proc );

#endif
