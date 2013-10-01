#include "config.h"

#include "auxiliary.h"
#include "iterators.h"

#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h" 
#include "starformation_feedback.h"
#include "step.h"
#include "hydro.h"


void setup_cell_feedback(int level){}

#ifdef HYDRO
void cell_feedback(int level, int time_multiplier) {
	int i;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	double t_next, dt_next;
	
	start_time( WORK_TIMER );
	
	setup_cell_feedback(level);
	dt_next = time_multiplier*dtl[level];
	t_next = tl[level] + dt_next;
	
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#ifdef STAR_FORMATION
#pragma omp parallel for default(none), private(iter_cell), shared(num_level_cells,level_cells,level,t_next, sf_feedback_cell), schedule(dynamic)
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];
		sf_feedback_cell->hydro_feedback_cell(level,iter_cell,t_next, dt_next);  
	}
	update_buffer_level( level, all_hydro_vars, num_hydro_vars );
#pragma omp parallel for default(none), private(iter_cell), shared(num_level_cells,level_cells,level,t_next, sf_feedback_cell), schedule(dynamic)
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];
		sf_feedback_cell->nonlocal_feedback_cell(level,iter_cell,t_next, dt_next);  
	}
#endif /* STAR_FORMATION*/
	
	cart_free(level_cells);
	
	end_time( WORK_TIMER );
}

#endif /* HYDRO */
