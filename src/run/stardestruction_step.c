#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rand.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"

#include "starformation_step.h"
#include "step.h"

#ifdef LOG_STAR_CREATION
#include "logging.h"
#endif

#ifdef HYDRO

void star_destruction(int level) {
  int i, j;
  int icell, idelete, ipart, ipart_prev;
  int num_level_cells;
  int *level_cells;
  double dt_eff;

  if(sf_feedback_particle->destroy_star_particle == NULL) return;

  start_time( WORK_TIMER );

#pragma omp parallel for default(none), private(icell,ipart,ipart_prev,idelete), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_id,star_particle_type,particle_species_indices,num_particle_species,particle_list_next, particle_list_prev, sf_feedback_particle), schedule(dynamic)
    for ( i = 0; i < num_level_cells; i++ ) {
	icell = level_cells[i];
	
	ipart = cell_particle_list[icell];
        while ( ipart != NULL_PARTICLE ) {
		ipart_prev = particle_list_prev[ipart] ;
		if ( particle_is_star(ipart) ){
		    //stellar_destruction(level,cell,ipart,&idelete );
                    idelete = sf_feedback_particle->destroy_star_particle(level,icell,ipart);
                    if(idelete == -1){ 
                        delete_particle(icell,ipart);
                        particle_free(ipart);
                        /* deleted current particle so go back one */
                        ipart = ipart_prev;
                    } 
		}

		/* go to next particle in list */
		if( ipart != NULL_PARTICLE ){  ipart = particle_list_next[ipart]; }
	    }
    }
    cart_free(level_cells);

  end_time( WORK_TIMER );
}

#endif /* HYDRO */
#endif /* STAR_FORMATION */
