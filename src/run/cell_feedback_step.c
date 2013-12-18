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
	int num_level_cells;
	int *level_cells;
	double t_next, dt_next;

#ifdef STAR_FORMATION	
	if ( sf_feedback_cell->hydro_feedback_cell == NULL && sf_feedback_cell->gather_feedback_cell == NULL ) return;
	
	start_time( WORK_TIMER );
	
	setup_cell_feedback(level);
	dt_next = time_multiplier*dtl[level];
	t_next = tl[level] + dt_next;
	
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
	
	if( sf_feedback_cell->hydro_feedback_cell != NULL ){
#pragma omp parallel for default(none), shared(num_level_cells,level_cells,level,t_next,dt_next, sf_feedback_cell), schedule(dynamic)
		for ( i = 0; i < num_level_cells; i++ ) {
			sf_feedback_cell->hydro_feedback_cell(level, level_cells[i], t_next, dt_next);  
		}
		cart_free(level_cells);
	}

	if( sf_feedback_cell->gather_feedback_cell != NULL ){
		/* update is only necessary for the gather operation if that form of feedback happens */
		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
#pragma omp parallel for default(none), shared(num_level_cells,level_cells,level,t_next,dt_next, sf_feedback_cell), schedule(dynamic)
		for ( i = 0; i < num_level_cells; i++ ) {
			sf_feedback_cell->gather_feedback_cell(level, level_cells[i], t_next, dt_next);  
		}
	}

	cart_free(level_cells);
	end_time( WORK_TIMER );
#endif /* STAR_FORMATION */

}

#endif /* HYDRO */
