#include "config.h"
#ifdef STARFORM

#include "auxiliary.h"
#include "iterators.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "step.h"


#ifdef BLASTWAVE_FEEDBACK
typedef struct
{
  double energy;
  double metals;
  double dt;
}
fb_pars;

extern fb_pars fbp_bw_phys;
#endif/* BLASTWAVE_FEEDBACK */


#ifdef BLASTWAVE_FEEDBACK
void check_bwtime_precision(int level)
{
  /* unlikely but tl and dtl are double */
  cart_assert( (dtl[level]*units->time/constants->yr) / fbp_bw_phys.dt > 1e-6 ); 
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
#ifdef STAR_PARTICLE_TYPES
#pragma omp parallel for default(none), private(iter_cell,ipart), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,t_next,particle_id,star_particle_type,particle_species_indices,num_particle_species,particle_list_next), schedule(dynamic)
#else
#pragma omp parallel for default(none), private(iter_cell,ipart), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,t_next,particle_id,particle_species_indices,num_particle_species,particle_list_next), schedule(dynamic)
#endif
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];

		ipart = cell_particle_list[iter_cell];
		while ( ipart != NULL_PARTICLE ) {
			if ( particle_is_star(ipart) && particle_t[ipart] < t_next - 0.5*dtl[max_level] ) {
#ifdef STAR_PARTICLE_TYPES
				if ( star_particle_type[ipart] == STAR_TYPE_NORMAL ) {                                                                              
					stellar_feedback(level,iter_cell,ipart,t_next);
				}
#else
				stellar_feedback(level,iter_cell,ipart,t_next);
#endif /* STAR_PARTICLE_TYPES */
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

#endif /* STARFORM */
