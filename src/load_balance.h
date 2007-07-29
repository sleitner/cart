#ifndef __LOAD_BALANCE_H__
#define __LOAD_BALANCE_H__

extern float cost_per_cell;
extern float cost_per_particle;
extern float est_buffer_fraction;
extern int load_balance_frequency;

int divide_list_recursive( float *global_work, int *global_counts, int num_root_cells_in_division,
	double total_work, long total_counts, int num_procs_in_division, int first_proc,
	int first_cell_index, int *proc_index );
int divide_list_linear( float *global_work, int *global_counts, int num_root_cells_in_division,
        double total_work, long total_counts, int num_procs_in_division, int first_proc,
        int first_cell_index, int *proc_index );
void load_balance_entire_volume( float *global_work, int *global_cells, int *new_proc_sfc_index );
void load_balance();

#endif
