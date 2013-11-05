#include "config.h"

#include "auxiliary.h"
#include "iterators.h"

#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h" 
#include "starformation_feedback.h"
#include "cell_buffer.h"
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
	
#ifdef STAR_FORMATION
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(iter_cell), shared(num_level_cells,level_cells,level,t_next,dt_next, sf_feedback_cell), schedule(dynamic)
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];
		sf_feedback_cell->hydro_feedback_cell(level,iter_cell,t_next, dt_next);  
	}
	cart_free(level_cells);
        /* update is only necessary for the gather operation if that form of feedback happens */
        update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(iter_cell), shared(num_level_cells,level_cells,level,t_next,dt_next, sf_feedback_cell), schedule(dynamic)
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];
		sf_feedback_cell->gather_feedback_cell(level,iter_cell,t_next, dt_next);  
	}
	cart_free(level_cells);
#endif /* STAR_FORMATION*/
	
	end_time( WORK_TIMER );
}

#endif /* HYDRO */
