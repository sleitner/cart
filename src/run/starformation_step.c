#include "config.h"
#if defined(STAR_FORMATION) && defined(HYDRO)

#include <math.h>

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rand.h"
#include "starformation.h"
#include "starformation_recipe.h"
#include "starformation_formstar.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "starformation_step.h"
#include "step.h"

#ifdef LOG_STAR_CREATION
#include "logging.h"
#endif

/* star formation parameters */
extern int sf_min_level;

void star_formation(int level, int time_multiplier) {
  int i, j;
  int icell;
  int num_level_cells;
  int *level_cells;
  double dt_eff;
  float *sfr;

  if ( level < sf_min_level ) return;

  start_time( WORK_TIMER );

  sfr = cart_alloc(float,num_level_cells);
  star_formation_rate(level,num_level_cells,level_cells,sfr);

  dt_eff = dtl[level] * time_multiplier;

  select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,level,sf_recipe,dt_eff, sfr, dtl), schedule(dynamic)
  for(i=0;i<num_level_cells; i++)
    {
        if ( sfr[i] <=0 ) continue;
        icell = level_cells[i];
        sf_formstar->form_star_particles(level,icell,dtl[level],dt_eff,sfr[i]);
    }

  cart_free(level_cells);
  cart_free(sfr);
  end_time( WORK_TIMER );
}

void remap_star_ids() {
	int i;
	int proc;
	int ipart;
	int block;
	int max_stars;	
	int new_id;
	int total_new_stars;
	int *block_ids;
	int proc_new_stars[MAX_PROCS];

	start_time( COMMUNICATION_TIMER );

	/* collect number of stars created */
	MPI_Allgather( &num_new_stars, 1, MPI_INT, proc_new_stars, 1, MPI_INT, mpi.comm.run );

	/* find how many "blocks" to expect */
	max_stars = 0;
	total_new_stars = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc_new_stars[proc] > max_stars ) {
			max_stars = proc_new_stars[proc];
		}

		total_new_stars += proc_new_stars[proc];
	}

	if ( total_new_stars > 0 ) {
		/* create lists of indices for each block */
		block_ids = cart_alloc(int, max_stars );

		block_ids[0] = 0;
		for ( block = 1; block < max_stars; block++ ) {
			block_ids[block] = block_ids[block-1];
			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( proc_new_stars[proc] >= block ) {
					block_ids[block]++;
				}
			}
		}
	
		/* find all newly allocated stars and remap their id's (keeping order) */
		for ( ipart = 0; ipart < num_star_particles; ipart++ ) {
			if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && 
					particle_id[ipart] >= particle_species_indices[num_particle_species] ) {
	
				block = ( particle_id[ipart] - particle_species_indices[num_particle_species] ) / num_procs;
				proc = ( particle_id[ipart] - particle_species_indices[num_particle_species] ) % num_procs;
				new_id = particle_species_indices[num_particle_species] + block_ids[block];
	
				for ( i = 0; i < proc; i++ ) {
					if ( proc_new_stars[i] > block ) {
							new_id++;
					}
				}
				
				cart_assert( new_id <= particle_id[ipart] && 
					new_id < particle_species_indices[num_particle_species]+total_new_stars );

				particle_id[ipart] = new_id;
			}
		}

		cart_free(block_ids);
	}

	particle_species_indices[num_particle_species] += total_new_stars;
	particle_species_num[num_particle_species-1] += total_new_stars;
	num_particles_total += total_new_stars;
	num_new_stars = 0;

	end_time( COMMUNICATION_TIMER );
}

#endif /* STAR_FORMATION && HYDRO */
