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
void cell_feedback(int level) {
	int i;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	double t_next;

	start_time( WORK_TIMER );

	setup_cell_feedback(level);
	t_next = tl[level] + dtl[level];

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#ifdef STAR_FORMATION
#pragma omp parallel for default(none), private(iter_cell), shared(num_level_cells,level_cells,level,t_next), schedule(dynamic)
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];
                sf_feedback_cell->hydro_feedback_cell(level,iter_cell,t_next, dtl[level]);  
	}
#endif /* STAR_FORMATION*/

	cart_free(level_cells);

	end_time( WORK_TIMER );
}

#endif /* HYDRO */
