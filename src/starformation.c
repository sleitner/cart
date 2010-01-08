#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_recipes.h"
#include "starformation_feedback.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"


int num_local_star_particles;
int last_star_id;
int num_new_stars;

double total_stellar_mass = 0.0;
double total_stellar_initial_mass = 0.0;

float star_tbirth[num_star_particles];
float star_initial_mass[num_star_particles];

#ifdef ENRICH
float star_metallicity_II[num_star_particles];
#ifdef ENRICH_SNIa
float star_metallicity_Ia[num_star_particles];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

float star_formation_volume_min[nDim];
float star_formation_volume_max[nDim];

DEFINE_LEVEL_ARRAY(int,star_formation_frequency);

/* star formation parameters */
int sf_recipe = 0;             /* default to the old-style recipe */
int sf_min_level = min_level;  /* Minimum level on which to create stars */

double sf_min_gas_number_density = 0.1;      /* in cm^{-3}; used to be called rho_SF */
double sf_max_gas_temperature = 2.0e4;       /* in K; used to be called T_SF */
double sf_timescale = 3.0e7;                 /* in yrs; used to be called tau_SF, did not exist in HART */
double sf_sampling_timescale = 1.0e6;        /* in yrs; used to be called dtmin_SF, also in HART */
double sf_min_stellar_particle_mass = 0.0;   /* in Msun; used to be called dm_star_min */


void config_init_star_formation()
{
  /*
  //  General parameters
  */
  control_parameter_add2(control_parameter_int,&sf_recipe,"sf:recipe","sf_recipe","recipe for star formation. Available recipes: \n   0 (oldstyle HART recipe),\n   1 (Gnedin et al 2009 recipes).\nUse <sf:recipe=1:min-cloud-density> and <sf:recipe=1:max-cloud-density> to mimic Gnedin et al 2009 recipes 1 to 3.");

  control_parameter_add2(control_parameter_int,&sf_min_level,"sf:min-level","sf_min_level","minimum level on which do star formation. Cells with level < <sf:min-level> form no stars no matter what.");

  control_parameter_add3(control_parameter_double,&sf_min_gas_number_density,"sf:min-gas-number-density","sf_min_gas_number_density","rho_sf","the gas total hydrogen number density threshold for star formation, in cm^{-3}. No star formation is done for gas at lower densities.");

  control_parameter_add3(control_parameter_double,&sf_max_gas_temperature,"sf:max-gas-temperature","sf_max_gas_temperature","t_sf","the maximum gas temperature (in K) for star formation. No star formation is done in hotter gas.");

  control_parameter_add3(control_parameter_double,&sf_timescale,"sf:timescale","sf_timescale","tau_sf","the timescale for star formation. Star formation in a given cell is assumed to continue with the constant rate for that period of time.");

  control_parameter_add3(control_parameter_double,&sf_sampling_timescale,"sf:sampling-timescale","sf_sampling_timescale","dtmin_sf","the timescale on which the conditions for star formation are checked. This is a numerical parameter only, no physical results should depend on it; its value should be sufficiently smaller than the <sf:timescale> parameter.  This parameter used to be called 'dtmin_SF' in HART.");

  control_parameter_add3(control_parameter_double,&sf_min_stellar_particle_mass,"sf:min-stellar-particle-mass","sf_min_stellar_particle_mass","dm_star_min","minimum mass for a newly created stellar particle, in solar masses. This value should be small enough to avoid artifically boosting the SFR in the low density gas.");

  config_init_star_formation_recipes();
  config_init_star_formation_feedback();
}


void config_verify_star_formation()
{
  /*
  //  General parameters
  */
  cart_assert(sf_recipe>=0 && sf_recipe<=1);

  cart_assert(sf_min_level>=min_level && sf_min_level<=max_level);

  cart_assert(sf_min_gas_number_density > 0.0);

  cart_assert(sf_max_gas_temperature > 10.0);

  cart_assert(sf_timescale > 0.0);

  cart_assert(sf_sampling_timescale < 0.5*sf_timescale);

  cart_assert(!(sf_min_stellar_particle_mass < 0.0));

  config_verify_star_formation_recipes();
  config_verify_star_formation_feedback();
}


void init_star_formation()
{
  int j;
  
  init_star_formation_feedback();

  /*
  //  NG: This limit is very obscure, disable by default
  */
  for(j=0; j<nDim; j++)
    {
      star_formation_volume_min[j] = 0.0;
      star_formation_volume_max[j] = num_grid;
    }
}


#ifdef HYDRO

void star_formation_rate(int level, int num_level_cells, int *level_cells, float *sfr)
{
  int i, j;
  int cell;
  float pos[nDim];
  int do_star_formation;
  double tem_max, rho_min;

  tem_max = sf_max_gas_temperature/(constants->wmu*units->temperature);
  rho_min = sf_min_gas_number_density/(constants->XH*units->number_density);
#ifdef COSMOLOGY
  rho_min = max(rho_min,200*cosmology->OmegaB/cosmology->OmegaM);
#endif

  setup_star_formation_recipes(level);

  for(i=0; i<num_level_cells; i++)
    {
      cell = level_cells[i];

      sfr[i] = -1.0;
      if(cell_is_leaf(cell))
	{
	  /* check position */
	  cell_position(cell,pos);

	  do_star_formation = 1;
	  for(j=0; j<nDim; j++)
	    {
	      if(pos[j]<star_formation_volume_min[j] || pos[j]>star_formation_volume_max[j])
		{
		  do_star_formation = 0;
		}
	    }

	  if(do_star_formation && cell_gas_density(cell)>rho_min && cell_gas_pressure(cell)/cell_gas_density(cell)<tem_max)
	    {
	      sfr[i] = sf_rate(cell);
	    }
	}
    }
}


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
  double P_SF;
  float *sfr;

  if ( level < sf_min_level) return;

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

      /* randomly generate particle on timescale tau_SF_eff */
      if ( cart_rand() > P_SF )
	{
	  dm_star = min( max(sfr[i]*dt_SF*cell_volume[level],dm_star_min), cell_fraction * cell_gas_density(icell) );

	  /* create the new star */
	  create_star_particle( icell, dm_star );
	}
    }

  cart_free( sfr );
  cart_free( level_cells );
}

/*
//  NG: I haven't modified Doug's code below.
*/
void create_star_particle( int icell, float mass ) {
	int i;
	int ipart;
	int id;
	int level;
	float pos[nDim];
	float new_density;
	float density_fraction;

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( mass > 0.0 );

	id = last_star_id + local_proc_id + 1;
	last_star_id += num_procs;
	num_new_stars++;

	ipart = particle_alloc( id );
	cart_assert( ipart < num_star_particles );

	/* place particle at center of cell with cell momentum */
	cell_position(icell, pos );
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

	cell_gas_density(icell) = new_density;
	cell_gas_energy(icell) *= density_fraction;
	cell_gas_internal_energy(icell) *= density_fraction;
	cell_gas_pressure(icell) *= density_fraction;
	cell_momentum(icell,0) *= density_fraction;
	cell_momentum(icell,1) *= density_fraction;
	cell_momentum(icell,2) *= density_fraction;
		
#ifdef ENRICH
	cell_gas_metal_density_II(icell) = max( 1e-17,
		cell_gas_metal_density_II(icell) - star_metallicity_II[ipart]*mass*cell_volume_inverse[level] );
#ifdef ENRICH_SNIa
	cell_gas_metal_density_Ia(icell) = max( 1e-17,
		cell_gas_metal_density_Ia(icell) - star_metallicity_Ia[ipart]*mass*cell_volume_inverse[level] );
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
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

	/* collect number of stars created */
	MPI_Allgather( &num_new_stars, 1, MPI_INT, proc_new_stars, 1, MPI_INT, MPI_COMM_WORLD );

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
}

#endif /* STARFORM */
