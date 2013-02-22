#include "config.h"
#ifdef STAR_FORMATION

#include "auxiliary.h"
#include "iterators.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "agn_step.h"
#include "starformation_feedback_step.h"
#include "step.h"

#ifdef BLASTWAVE_FEEDBACK
extern double blastwave_time;
#endif/* BLASTWAVE_FEEDBACK */


#ifdef BLASTWAVE_FEEDBACK
void check_bwtime_precision(int level)
{
  /* unlikely but tl and dtl are double */
  cart_assert( (dtl[level]*units->time/constants->yr) / blastwave_time > 1e-6 ); 
}
#endif /* BLASTWAVE_FEEDBACK */

void star_particle_feedback(int level) {
	int i;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	double t_next;

	start_time( WORK_TIMER );

	setup_star_formation_feedback(level);
	t_next = tl[level] + dtl[level];

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#ifndef COMPILER_GCC
		/* Get compiler segfault under GCC */
#ifdef STAR_PARTICLE_TYPES
#pragma omp parallel for default(none), private(iter_cell,ipart), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,t_next,particle_id,star_particle_type,particle_species_indices,num_particle_species,particle_list_next, particle_mass), schedule(dynamic) 
#else
#pragma omp parallel for default(none), private(iter_cell,ipart), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,t_next,particle_id,particle_species_indices,num_particle_species,particle_list_next), schedule(dynamic)
#endif
#endif
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];

		ipart = cell_particle_list[iter_cell];
		while ( ipart != NULL_PARTICLE ) {
			if ( particle_is_star(ipart) && particle_t[ipart] < t_next - 0.5*dtl[max_level] 
#ifdef STAR_PARTICLE_TYPES
			     && (   star_particle_type[ipart] == STAR_TYPE_NORMAL 
                                 || star_particle_type[ipart] == STAR_TYPE_STARII 
                                 || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH) 
#endif /* STAR_PARTICLE_TYPES */
			) {
				stellar_feedback(level,iter_cell,ipart,t_next);
			}

			ipart = particle_list_next[ipart];
		}
	}

	cart_free(level_cells);

#if defined(STAR_PARTICLE_TYPES) && defined(AGN)
	agn_feedback( level );
#endif /* STAR_PARTICLE_TYPES && AGN */

	end_time( WORK_TIMER );
}

#endif /* STAR_FORMATION */
