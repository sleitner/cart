#ifndef __LOAD_BALANCE_H__
#define __LOAD_BALANCE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


extern float cost_per_cell;
extern float cost_per_particle;
extern float est_buffer_fraction;
extern int load_balance_frequency;

#ifdef PARTICLES
#define num_constraints         2
#else
#define num_constraints         1
#endif

void config_init_load_balance();
void config_verify_load_balance();

int divide_list_recursive( float *global_work,
                int *constrained_quantities,
                int *per_proc_constraints,
                int num_root_cells_in_division, double total_work,
                int num_procs_in_division, int first_proc,
                int first_cell_index, int *proc_index );
int divide_list_linear( float *global_work, int *constrained_quantities,
                int *per_proc_constraints,
                int num_root_cells_in_division, double total_work,
                int num_procs_in_division, int first_proc,
                int first_cell_index, int *proc_index );
void load_balance_entire_volume( float *global_work, 
		int *constrained_quantities, 
		int *new_proc_sfc_index );
void load_balance();

#endif
