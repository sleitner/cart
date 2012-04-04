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
#include "units.h"

#include "starformation_step.h"
#include "step.h"

#ifdef LOG_STAR_CREATION
#include "logging.h"
#endif


/* star formation parameters */
extern int sf_min_level;
extern float sf_metallicity_floor;
extern double sf_timescale;
extern double sf_min_stellar_particle_mass;


#ifdef HYDRO

void star_formation( int level, int time_multiplier )
{
  int i;
  int icell;
  int num_level_cells;
  int *level_cells;
  double dt_SF;
  double dt_eff, mstar_min;
  double mstar;
#ifdef OLDSTYLE_SF_ALGORITHM
  double P_SF, P_mass;
#endif
  float *sfr;

  if ( level < sf_min_level ) return;

  start_time( WORK_TIMER );

  mstar_min = sf_min_stellar_particle_mass * constants->Msun / units->mass; 

  dt_SF = sf_timescale * constants->yr / units->time;
  dt_eff = dtl[level] * time_multiplier;

#ifdef OLDSTYLE_SF_ALGORITHM
  /* probability of forming a star is Poisson with <t> = dt_SF */
  P_SF = exp( -dt_eff / dt_SF );
#endif /* OLDSTYLE_SF_ALGORITHM */

  select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
  sfr = cart_alloc(float,num_level_cells);
  star_formation_rate(level,num_level_cells,level_cells,sfr);

  for ( i = 0; i < num_level_cells; i++ ) { 
	  if ( sfr[i] > 0.0 ) {
		icell = level_cells[i];
		mstar = sfr[i]*dt_SF*cell_volume[level];

#ifdef OLDSTYLE_SF_ALGORITHM

		if(mstar < mstar_min) P_mass = 1 - mstar/mstar_min; else P_mass = 0;

		/* randomly generate particle on timescale dt_SF */
		if(cart_rand() > P_SF+P_mass-P_SF*P_mass)
		  mstar = max(mstar_min,mstar);
		else 
		  mstar = 0.0;
		  
#else  /* OLDSTYLE_SF_ALGORITHM */

		/* draw number of star formation events 0...\inf from poisson distribution */
		mstar = max( mstar_min, mstar );
		mstar *= (double)cart_rand_poisson( sfr[i]*cell_volume[level]*dt_eff/mstar );

#endif /* OLDSTYLE_SF_ALGORITHM */

		if ( mstar > 0.0 ) {

			/* create the new star */
#ifdef LOG_STAR_CREATION	  
		  	log_star_creation( icell, mstar, FILE_RECORD);
#endif
	   		create_star_particle( icell, mstar, dtl[level], STAR_TYPE_NORMAL );

#ifdef BLASTWAVE_FEEDBACK
			init_blastwave(icell);
#endif /* BLASTWAVE_FEEDBACK */
		}
	 }
  }

  cart_free( sfr );
  cart_free( level_cells );

  end_time( WORK_TIMER );
}

#endif /* HYDRO */

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

#endif /* STAR_FORMATION */
