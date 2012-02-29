#ifndef __PACK_H__
#define __PACK_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "skiplist.h"

#define	PACK_COMPLETED		(-1)

typedef struct PACK {
	int cell_type;
	skiplist *tree_list[MAX_PROCS];

	/* arrays for sending */
	int num_sending_cells_total[MAX_PROCS];
	int num_sending_cells[MAX_PROCS][max_level-min_level+1];

	int *root_cells[MAX_PROCS];
	int *cell_refined[MAX_PROCS];
	float *cell_vars[MAX_PROCS];

	/* arrays for receiving */
	int num_receiving_cells_total[MAX_PROCS];
        int num_receiving_cells[MAX_PROCS][max_level-min_level+1];
} pack;

pack *pack_init( int cell_type );
void pack_destroy( pack *p );
void pack_add_root_trees( pack *p, int *new_proc_sfc_index, int sfc1, int sfc2 );
void pack_add_root_tree( pack *p, int proc, int sfc );
void pack_apply( pack *p );
void pack_communicate( pack *p );

#endif
