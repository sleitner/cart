#include "config.h"
#ifdef STARFORM

#include "agn.h"
#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "starformation.h"
#include "starformation_recipes.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"


int num_local_star_particles = 0;
int last_star_id = -1;
int num_new_stars = 0;

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
double sf_timescale = 3.0e7;                 /* in yrs; used to be called tau_SF, did not exist in HART */
double sf_sampling_timescale = 1.0e6;        /* in yrs; used to be called dtmin_SF, also in HART */
double sf_min_stellar_particle_mass = 0.0;   /* in Msun; used to be called dm_star_min */

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

  control_parameter_add3(control_parameter_double,&sf_timescale,"sf:timescale","sf_timescale","tau_SF","the timescale for star formation. Star formation in a given cell is assumed to continue with the constant rate for that period of time.");

  control_parameter_add3(control_parameter_double,&sf_sampling_timescale,"sf:sampling-timescale","sf_sampling_timescale","dtmin_sf","the timescale on which the conditions for star formation are checked. This is a numerical parameter only, no physical results should depend on it; its value should be sufficiently smaller than the <sf:timescale> parameter.  This parameter used to be called 'dtmin_SF' in HART.");

  control_parameter_add3(control_parameter_double,&sf_min_stellar_particle_mass,"sf:min-stellar-particle-mass","sf_min_stellar_particle_mass","dm_star_min","minimum mass for a newly created stellar particle, in solar masses. This value should be small enough to avoid artifically boosting the SFR in the low density gas.");

  control_parameter_add2(control_parameter_float,&sf_min_overdensity,"sf:min-overdensity","sf_min_overdensity","the value of the overdensity (in total mass, not gas mass) below which the gas density threshold is not allowed to descend. I.e., rho_min = max(<sf:min-gas-number-density>/(constants->XH*units->number_density),<sf:min-overdensity>*cosmology->OmegaB/cosmology->OmegaM");

  control_parameter_add2(control_parameter_float,&sf_metallicity_floor,"sf:metallicity-floor","sf_metallicity_floor","the minimum amount of metallicity (in solar units) sprinkled in the cell where a stellar particle is created. This parameter is designed to emulate metal enrichment from the first stars and/or early, unresolved episode of star formation. Just BEFORE a stellar particle is created, if the metallicity in the cell is less than this value, it is increased to this value. Thus, all stellar particles are created in the gas with at least this metallicity. This value is NOT taken into account when computing the SFR, only when actually creating a stellar particle. It is not related to rt:dust-to-gas floor, which MAY affect the SFR.");

  config_init_star_formation_recipes();
  config_init_star_formation_feedback();
}


void config_verify_star_formation()
{
  /*
  //  General parameters
  */
  cart_assert(sf_min_level>=min_level && sf_min_level<=max_level);

  cart_assert(sf_min_gas_number_density > 0.0);

  cart_assert(sf_max_gas_temperature > 10.0);

  cart_assert(sf_timescale > 0.0);

  cart_assert(sf_sampling_timescale < 0.5*sf_timescale);

  cart_assert(!(sf_min_stellar_particle_mass < 0.0));

  cart_assert(!(sf_metallicity_floor < 0.0));


  config_verify_star_formation_recipes();
  config_verify_star_formation_feedback();
}

#endif /* HYDRO */


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
  double pos[nDim];
  int do_star_formation;
  double tem_max, rho_min;

  tem_max = sf_max_gas_temperature/(constants->wmu*units->temperature);
  rho_min = sf_min_gas_number_density/(constants->XH*units->number_density);
#ifdef COSMOLOGY
  rho_min = max(rho_min,sf_min_overdensity*cosmology->OmegaB/cosmology->OmegaM);
#endif

  sf_recipe->setup(level);

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

#endif /* HYDRO */

#endif /* STARFORM */
