#include "config.h"
#ifdef STARFORM

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


void create_star_particle( int icell, float mass, int type );


#ifdef HYDRO

void star_formation( int level, int time_multiplier )
{
  int i;
  int icell;
  int num_level_cells;
  int *level_cells;
  double cell_fraction;
  double dm_star;
  double dt_SF;
  double dt_eff, dm_star_min;
  double P_SF, P_mass, mstar;
  float *sfr;

  if ( level < sf_min_level) return;

  start_time( WORK_TIMER );

  cell_fraction = 0.667 * cell_volume[level];
  dm_star_min = sf_min_stellar_particle_mass * constants->Msun / units->mass; 

  dt_SF = sf_timescale * constants->yr / units->time;
  dt_eff = dtl[level] * time_multiplier;

  /* probability of forming a star is Poisson with <t> = tau_SF */
  P_SF = exp( -dt_eff / dt_SF );

  select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
  sfr = cart_alloc(float,num_level_cells);

  star_formation_rate(level,num_level_cells,level_cells,sfr);

  for ( i = 0; i < num_level_cells; i++ ) if ( sfr[i] > 0.0 )
    {
      icell = level_cells[i];

      /*
      // NEW PART ADDED BY OLEG
      */

      mstar = sfr[i]*dt_SF*cell_volume[level];
      if(mstar < dm_star_min) P_mass = 1 - mstar/dm_star_min; else P_mass = 0;

      /* randomly generate particle on timescale tau_SF_eff */
      //      if(cart_rand()>P_SF && cart_rand()>P_mass)
      if(cart_rand() > P_SF+P_mass-P_SF*P_mass)
	{
	  dm_star = min( max(mstar,dm_star_min), cell_fraction * cell_gas_density(icell) );

	  /*
	  // NEW PART ENDS HERE
	  */

	  /* create the new star */
#ifdef LOG_STAR_CREATION	  
	  log_star_creation( icell, dm_star, FILE_RECORD);
#endif
	  create_star_particle( icell, dm_star, STAR_TYPE_NORMAL );

#ifdef BLASTWAVE_FEEDBACK
	  init_blastwave(icell);
#endif /* BLASTWAVE_FEEDBACK */
	}
    }

  cart_free( sfr );
  cart_free( level_cells );

  end_time( WORK_TIMER );
}

/*
//  NG: I haven't modified Doug's code below.
*/
void create_star_particle( int icell, float mass, int type ) {
	int i;
	int ipart;
	int id;
	int level;
	double pos[nDim];
	float new_density;
	float density_fraction, thermal_pressure;

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( mass > 0.0 );

	id = last_star_id + local_proc_id + 1;
	last_star_id += num_procs;
	num_new_stars++;

	ipart = particle_alloc( id );
	cart_assert( ipart < num_star_particles );

	/*
	//  This is an obscure parameter, read its help string in 
	//  config_init_star_formation().
	*/
#ifdef ENRICH
	if(sf_metallicity_floor>0.0 && cell_gas_metal_density_II(icell)<sf_metallicity_floor*constants->Zsun*cell_gas_density(icell))
	  {
	    cell_gas_metal_density_II(icell) =  sf_metallicity_floor*constants->Zsun*cell_gas_density(icell);
	  }
#endif

#ifdef STAR_PARTICLE_TYPES
	star_particle_type[ipart] = type;
#endif

	/* place particle at center of cell with cell momentum */
	cell_center_position(icell, pos );
	level = cell_level(icell);

	for ( i = 0; i < nDim; i++ ) {
		particle_x[ipart][i] = pos[i];
	}

	for ( i = 0; i < nDim; i++ ) {
		particle_v[ipart][i] = cell_momentum(icell,i) / cell_gas_density(icell);
	}

	particle_t[ipart] = tl[level];
	particle_dt[ipart] = dtl[level];

	star_tbirth[ipart] = tl[level];
	particle_mass[ipart] = mass;
	star_initial_mass[ipart] = mass;

#ifdef ENRICH
	star_metallicity_II[ipart] = cell_gas_metal_density_II(icell) / cell_gas_density(icell);
#ifdef ENRICH_SNIa
	star_metallicity_Ia[ipart] = cell_gas_metal_density_Ia(icell) / cell_gas_density(icell);
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
	
	/* insert particle into cell linked list */
	insert_particle( icell, ipart );

	/* adjust cell values */
	new_density = cell_gas_density(icell) - mass * cell_volume_inverse[level];
	density_fraction = new_density / cell_gas_density(icell);

	/*
	// NG: this is to allow non-thermal pressure contribution
	*/
	thermal_pressure = max((cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell),0.0);
	cell_gas_pressure(icell) = max(0.0,cell_gas_pressure(icell)-thermal_pressure);

	cell_gas_density(icell) = new_density;
	cell_gas_energy(icell) *= density_fraction;
	cell_gas_internal_energy(icell) *= density_fraction;
	cell_momentum(icell,0) *= density_fraction;
	cell_momentum(icell,1) *= density_fraction;
	cell_momentum(icell,2) *= density_fraction;
		
	cell_gas_pressure(icell) += thermal_pressure*density_fraction;

	for ( i = 0; i < num_chem_species; i++ ) {
		cell_advected_variable(icell,i) *= density_fraction;
	}
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

#endif /* STARFORM */
