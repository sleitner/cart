#include "config.h"
#if defined(PARTICLES) && defined(STAR_FORMATION)

#include "agn.h"
#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_recipe.h"
#include "starformation_feedback.h"
#include "starformation_formstar.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "imf.h"

int num_local_star_particles = 0;
particleid_t last_star_id = -1;
int num_new_stars = 0;

double total_stellar_mass = 0.0;
double total_stellar_initial_mass = 0.0;

float star_tbirth[num_star_particles];
float star_initial_mass[num_star_particles];

#ifdef ENRICHMENT
float star_metallicity_II[num_star_particles];
#ifdef ENRICHMENT_SNIa
float star_metallicity_Ia[num_star_particles];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

#ifdef STAR_PARTICLE_TYPES
int star_particle_type[num_star_particles];
#endif /* STAR_PARTICLE_TYPES */

float star_formation_volume_min[nDim];
float star_formation_volume_max[nDim];


#ifdef HYDRO

DEFINE_LEVEL_ARRAY(int,star_formation_frequency);

/* star formation parameters */
int sf_min_level = min_level;  /* Minimum level on which to create stars */

double sf_min_gas_number_density = 0.1;      /* in cm^{-3}; used to be called rho_SF */
double sf_max_gas_temperature = 2.0e4;       /* in K; used to be called T_SF */

double sf_sampling_timescale = 1.0e6;        /* in yrs; used to be called dtmin_SF, also in HART */

float sf_min_overdensity = 200;

float sf_metallicity_floor = 0.0;            /* this is an obscure parameter, read its help string in config_init_star_formation(). */

void config_init_star_formation()
{
  /*
  //  General parameters
  */
  control_parameter_add2(control_parameter_int,&sf_min_level,"sf:min-level","sf_min_level","minimum level on which do star formation. Cells with level < <sf:min-level> form no stars no matter what.");

  control_parameter_add3(control_parameter_double,&sf_min_gas_number_density,"sf:min-gas-number-density","sf_min_gas_number_density","rho_sf","the gas total hydrogen number density threshold for star formation, in cm^{-3}. No star formation is done for gas at lower densities.");

  control_parameter_add3(control_parameter_double,&sf_max_gas_temperature,"sf:max-gas-temperature","sf_max_gas_temperature","t_sf","the maximum gas temperature (in K) for star formation. No star formation is done in hotter gas.");

  control_parameter_add3(control_parameter_time,&sf_sampling_timescale,"sf:sampling-timescale","sf_sampling_timescale","dtmin_sf","the timescale on which the conditions for star formation are checked. This is a numerical parameter only, no physical results should depend on it; its value should be sufficiently smaller than the <sf:timescale> parameter.  This parameter used to be called 'dtmin_SF' in HART.");

  control_parameter_add2(control_parameter_float,&sf_min_overdensity,"sf:min-overdensity","sf_min_overdensity","the value of the overdensity (in total mass, not gas mass) below which the gas density threshold is not allowed to descend. I.e., rho_min = MAX(<sf:min-gas-number-density>/(constants->XH*units->number_density),<sf:min-overdensity>*cosmology->OmegaB/cosmology->OmegaM");

  control_parameter_add2(control_parameter_float,&sf_metallicity_floor,"sf:metallicity-floor","sf_metallicity_floor","the minimum amount of metallicity (in solar units) sprinkled in the cell where a stellar particle is created. This parameter is designed to emulate metal enrichment from the first stars and/or early, unresolved episode of star formation. Just BEFORE a stellar particle is created, if the metallicity in the cell is less than this value, it is increased to this value. Thus, all stellar particles are created in the gas with at least this metallicity. This value is NOT taken into account when computing the SFR, only when actually creating a stellar particle. It is not related to rt:dust-to-gas floor, which MAY affect the SFR.");

  /* 
  //  IMF 
  */
  config_init_imf();

  config_init_star_formation_recipe();
  config_init_formstar();
  config_init_star_formation_feedback();
}


void config_verify_star_formation()
{
  /*
  //  General parameters
  */
  VERIFY(sf:min-level, sf_min_level>=min_level && sf_min_level<=max_level );

  VERIFY(sf:min-gas-number-density, sf_min_gas_number_density > 0.0 );

  VERIFY(sf:max-gas-temperature, sf_max_gas_temperature > 10.0 );

  VERIFY(sf:sampling-timescale, sf_sampling_timescale >= 0.0 ); 

  VERIFY(sf:metallicity-floor, !(sf_metallicity_floor < 0.0) );

  VERIFY(sf:min-overdensity, 1 );

  /* 
  //  IMF 
  */
  config_verify_imf();

  config_verify_star_formation_recipe();
  config_verify_star_formation_feedback();
  config_verify_formstar();
}

#endif /* HYDRO */


void init_star_formation()
{
  int j;
  
  init_star_formation_feedback();
  init_formstar();

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
  double pos[nDim];
  int do_star_formation;
  double tem_max, rho_min;

  tem_max = sf_max_gas_temperature/(constants->wmu*units->temperature);
  rho_min = sf_min_gas_number_density/(constants->XH*units->number_density);
#ifdef COSMOLOGY
  rho_min = MAX(rho_min,sf_min_overdensity*cosmology->OmegaB/cosmology->OmegaM);
#endif

  if(sf_recipe->setup != NULL) sf_recipe->setup(level);

  for(i=0; i<num_level_cells; i++)
    {
      cell = level_cells[i];

      sfr[i] = -1.0;
      if(cell_is_leaf(cell))
	{
	  /* check position */
	  cell_center_position(cell,pos);

	  do_star_formation = 1;
	  for(j=0; j<nDim; j++)
	    {
	      if(pos[j]<star_formation_volume_min[j] || pos[j]>star_formation_volume_max[j])
		{
		  do_star_formation = 0;
		}
	    }

	  if(do_star_formation && cell_gas_density(cell)>rho_min && cell_gas_temperature(cell)<tem_max)
	    {
	      sfr[i] = sf_recipe->rate(cell);
	    }
	}
    }
}

int create_star_particle( int icell, float mass, double pdt, int type ) {
	int i;
	int ipart;
	particleid_t id;
	int level;
	double pos[nDim];
	float new_density;
	float density_fraction, thermal_pressure;
	float true_mass;

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( mass > 0.0 );

	id = last_star_id + local_proc_id + 1;
	last_star_id += num_procs;
	num_new_stars++;

	ipart = particle_alloc( id );
	cart_assert( ipart < num_star_particles );

	/* ensure star particle cannot consume entire cell gas mass */
	true_mass = MIN( mass, 0.667*cell_volume[cell_level(icell)]*cell_gas_density(icell) );

	/*
	//  This is an obscure parameter, read its help string in 
	//  config_init_star_formation().
	*/
#ifdef ENRICHMENT
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
	particle_dt[ipart] = pdt;

	star_tbirth[ipart] = tl[level];
	particle_mass[ipart] = true_mass;
	star_initial_mass[ipart] = true_mass;

#ifdef STAR_PARTICLE_TYPES
	if( star_particle_type[ipart] == STAR_TYPE_NORMAL || star_particle_type[ipart] == STAR_TYPE_STARII || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH ) {
#endif /* STAR_PARTICLE_TYPES */

#ifdef ENRICHMENT
            star_metallicity_II[ipart] = cell_gas_metal_density_II(icell) / cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
            star_metallicity_Ia[ipart] = cell_gas_metal_density_Ia(icell) / cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

#ifdef STAR_PARTICLE_TYPES
	} else if( star_particle_type[ipart] == STAR_TYPE_AGN ) {
#ifdef ENRICHMENT
            star_metallicity_II[ipart] = 0.0;
#ifdef ENRICHMENT_SNIa
            star_metallicity_Ia[ipart] = 0.0;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
	} else {
            cart_error("Initial metallicity for star_particle_type %d must be explicitly defined on creation",star_particle_type[ipart]);
	}	
#endif /* STAR_PARTICLE_TYPES */

	/* insert particle into cell linked list */
	insert_particle( icell, ipart );

	/* adjust cell values */
	new_density = cell_gas_density(icell) - true_mass * cell_volume_inverse[level];
	density_fraction = new_density / cell_gas_density(icell);

	/*
	// NG: this is to allow non-thermal pressure contribution
	*/
	thermal_pressure = MAX((cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell),0.0);
	cell_gas_pressure(icell) = MAX(0.0,cell_gas_pressure(icell)-thermal_pressure);

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

#ifdef BLASTWAVE_FEEDBACK
	start_blastwave(icell);
#endif /* BLASTWAVE_FEEDBACK */

	return ipart;
}

#endif /* HYDRO */
#endif /* PARTICLES && STAR_FORMATION */
