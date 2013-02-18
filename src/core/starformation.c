#include "config.h"
#if defined(PARTICLES) && defined(STAR_FORMATION)
#include <math.h>
#include <stdio.h>

#include "agn.h"
#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_recipe.h"
#include "starformation_feedback.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "imf.h"
#include "rand.h"
#include "onestarfits.h"

int num_local_star_particles = 0;
int last_star_id = -1;
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
double sf_timescale = 3.0e7;                 /* in yrs; used to be called tau_SF, did not exist in HART */
double sf_sampling_timescale = 1.0e6;        /* in yrs; used to be called dtmin_SF, also in HART */
double sf_min_stellar_particle_mass = 0.0;   /* in Msun; used to be called dm_star_min */

float sf_min_overdensity = 200;
float sf_metallicity_floor = 0.0;            /* this is an obscure parameter, read its help string in config_init_star_formation(). */

double fast_growth_probability=0;            /* how often does efficient star formation happen? */
double fast_growth_multiplier=10;            /* how quickly does the efficient star formation happen? */
int cluster_buildup_indicator=0;             /* are we growing clusters or forming them instantly? */
double cluster_age_spread=15.0e6;            /* in yrs; age above which "cluster" formation ceases */
double cluster_min_expected_mass=100;        /* do not form a cluster particle if it is expected to grow to less than this mass */

/* STARII related */
void get_com_pos( int icell, double com_pos[nDim], int level );
double starII_highmass_slope=-2.35;          /* IMF slope used for sampling starII masses */
double starII_minimum_mass=8.0;              /* in Msun; mass above which individual stars are sampled */
int starII_indicator=0 ;                     /* are we forming individual massive stars? */
int starII_runaway_indicator=0;              /* are we letting individual massive stars runaway? */

extern struct StellarFeedback sf_feedback_PopM;


void config_init_star_formation()
{
  /*
  //  General parameters
  */
  control_parameter_add2(control_parameter_int,&sf_min_level,"sf:min-level","sf_min_level","minimum level on which do star formation. Cells with level < <sf:min-level> form no stars no matter what.");

  control_parameter_add3(control_parameter_double,&sf_min_gas_number_density,"sf:min-gas-number-density","sf_min_gas_number_density","rho_sf","the gas total hydrogen number density threshold for star formation, in cm^{-3}. No star formation is done for gas at lower densities.");

  control_parameter_add3(control_parameter_double,&sf_max_gas_temperature,"sf:max-gas-temperature","sf_max_gas_temperature","t_sf","the maximum gas temperature (in K) for star formation. No star formation is done in hotter gas.");

  control_parameter_add3(control_parameter_time,&sf_timescale,"sf:timescale","sf_timescale","tau_SF","the timescale for star formation. Star formation in a given cell is assumed to continue with the constant rate for that period of time.");

  control_parameter_add3(control_parameter_time,&sf_sampling_timescale,"sf:sampling-timescale","sf_sampling_timescale","dtmin_sf","the timescale on which the conditions for star formation are checked. This is a numerical parameter only, no physical results should depend on it; its value should be sufficiently smaller than the <sf:timescale> parameter.  This parameter used to be called 'dtmin_SF' in HART.");

  control_parameter_add3(control_parameter_double,&sf_min_stellar_particle_mass,"sf:min-stellar-particle-mass","sf_min_stellar_particle_mass","dm_star_min","minimum mass for a newly created stellar particle, in solar masses. This value should be small enough to avoid artifically boosting the SFR in the low density gas.");

  control_parameter_add2(control_parameter_float,&sf_min_overdensity,"sf:min-overdensity","sf_min_overdensity","the value of the overdensity (in total mass, not gas mass) below which the gas density threshold is not allowed to descend. I.e., rho_min = MAX(<sf:min-gas-number-density>/(constants->XH*units->number_density),<sf:min-overdensity>*cosmology->OmegaB/cosmology->OmegaM");

  control_parameter_add2(control_parameter_float,&sf_metallicity_floor,"sf:metallicity-floor","sf_metallicity_floor","the minimum amount of metallicity (in solar units) sprinkled in the cell where a stellar particle is created. This parameter is designed to emulate metal enrichment from the first stars and/or early, unresolved episode of star formation. Just BEFORE a stellar particle is created, if the metallicity in the cell is less than this value, it is increased to this value. Thus, all stellar particles are created in the gas with at least this metallicity. This value is NOT taken into account when computing the SFR, only when actually creating a stellar particle. It is not related to rt:dust-to-gas floor, which MAY affect the SFR.");

  /* cluster buildup parameters*/ 
  control_parameter_add2(control_parameter_bool, &cluster_buildup_indicator, "cluster:buildup-indicator", "cluster_buildup_indicator", "turn on cluster buildup module for star formation");

  control_parameter_add3(control_parameter_time, &cluster_age_spread,"cluster:age-spread","cluster_age_spread","cluster_age_spread","timescale over which star particles representing clusters are allowed to grow.");

  control_parameter_add2(control_parameter_double, &cluster_min_expected_mass,"cluster:min-expected-mass","cluster_min_expected_mass","the minimum mass expected from sfr*cluster_age_spread allowed to seed a cluster.");

  control_parameter_add2(control_parameter_double, &fast_growth_probability, "cluster:fast_growth_probability", "fast_growth_probability", "probability of a fast growth particle");
  control_parameter_add2(control_parameter_double, &fast_growth_multiplier, "cluster:fast_growth_multiplier", "fast_growth_multiplier", "multiplier for fast growth particle star formation rates");

  /* 
  /* STARII parameters*/ 
  control_parameter_add2(control_parameter_bool, &starII_indicator, "starII:indicator", "starII_indicator", "turn on starII star-formation");

  control_parameter_add2(control_parameter_bool, &starII_runaway_indicator, "starII:runaway-indicator", "starII_runaway_indicator", "turn on runaway starIIs");

  control_parameter_add2(control_parameter_double,&starII_highmass_slope,"starII:highmass-slope","starII_highmass_slope","IMF slope used for sampling starII masses.");

  control_parameter_add2(control_parameter_double,&starII_minimum_mass,"starII:minimum-mass","starII_minimum_mass","the minimum mass of 'virtual' starII particles in Msun.");

  /* 
  //  IMF 
  */
  config_init_imf();

  config_init_star_formation_recipe();
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

  VERIFY(sf:timescale, sf_timescale >= 0.0 );

  VERIFY(sf:sampling-timescale, sf_sampling_timescale <= 0.5*sf_timescale );

  VERIFY(sf:min-stellar-particle-mass, !(sf_min_stellar_particle_mass < 0.0) );

  VERIFY(sf:metallicity-floor, !(sf_metallicity_floor < 0.0) );

  VERIFY(sf:min-overdensity, 1 );

  VERIFY(cluster:fast_growth_multiplier, fast_growth_multiplier >= 0);
  VERIFY(cluster:fast_growth_probability, fast_growth_probability >= 0 && fast_growth_probability <=1 );

/* STARII */
  if(starII_indicator){
#ifndef STAR_PARTICLE_TYPES
      cart_error("STAR_PARTICLE_TYPES must be defined for starII_indicator True");
#endif /* STAR_PARTICLE_TYPES */
      VERIFY(cluster:buildup-indicator,cluster_buildup_indicator==1); /* starII requires cluster build up currently */
      VERIFY(starII:runaway-indicator,starII_runaway_indicator==1 || starII_runaway_indicator==0);
  }else{
      VERIFY(starII:runaway-indicator, starII_runaway_indicator==0);
  }
  VERIFY(starII:highmass-slope,starII_highmass_slope);
  VERIFY(starII:minimum-mass, starII_minimum_mass >1.0 );
  if(cluster_buildup_indicator){
      VERIFY(sf:min-stellar-particle-mass, sf_min_stellar_particle_mass == 0.0 ); 
      VERIFY(sf:timescale, sf_timescale == 0.0 );
      VERIFY(sf:sampling-timescale, sf_sampling_timescale == 0 );
  }
  VERIFY(cluster:age-spread, cluster_age_spread >1.0e6 );
  VERIFY(cluster:min-expected-mass, cluster_min_expected_mass >= 0.0 );
 
  /* 
  //  IMF 
  */
  config_verify_imf();

  config_verify_star_formation_recipe();
  config_verify_star_formation_feedback();
}

#endif /* HYDRO */

double mfrac_starII, starII_avg_mass, m_imf_tot;
double sIIminpw, imfmaxpw;
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

  double lowmass_avg_mass;
  double Nfrac_lowmass, Nfrac_starII;
  double N_imf_tot, N_imf_lowmass, N_imf_starII ;
  double m_imf_lowmass, m_imf_starII; 
  double mfrac_lowmass;
  
  if(starII_indicator){
      N_imf_tot = integrate( imf->f, imf->min_mass, imf->max_mass, 1e-6, 1e-9 );
      N_imf_lowmass = integrate( imf->f, imf->min_mass,starII_minimum_mass , 1e-6, 1e-9 );
      N_imf_starII = integrate( imf->f, starII_minimum_mass , imf->max_mass,1e-6, 1e-9 ); 
      Nfrac_lowmass = N_imf_lowmass/N_imf_tot;
      Nfrac_starII = (N_imf_tot - N_imf_lowmass)/N_imf_tot;
      
      m_imf_lowmass = integrate( imf->fm, imf->min_mass,starII_minimum_mass , 1e-6, 1e-9 );
      m_imf_starII = integrate( imf->fm, starII_minimum_mass, imf->max_mass , 1e-6, 1e-9 );
      m_imf_tot = integrate( imf->fm, imf->min_mass,imf->max_mass , 1e-6, 1e-9 );
      
      
      /*     mfrac_starII = m_imf_lowmass/N_imf_lowmass; */
      mfrac_starII = m_imf_starII/m_imf_tot;
      mfrac_lowmass = m_imf_lowmass/m_imf_tot;
      lowmass_avg_mass = m_imf_lowmass/N_imf_lowmass;
      starII_avg_mass = m_imf_starII/N_imf_starII;
      
      cart_debug("minimum %e Nfrac_lowmass %e Nfrac_starII %e ", 
		 starII_minimum_mass, Nfrac_lowmass, Nfrac_starII );
      cart_debug("mfrac_lowmass %e mfrac_starII %e ", 
	     mfrac_lowmass, mfrac_starII);
      cart_debug("avg lowmass star mass = %e ; avg starII mass %e", 
		 lowmass_avg_mass, starII_avg_mass);
      
      sIIminpw = pow( starII_minimum_mass, starII_highmass_slope+1 );
      imfmaxpw = pow( imf->max_mass, starII_highmass_slope+1 ) ;

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

double sample_exponential(double tau){
    return -tau*log(cart_rand());
}
double runaway_mean_kick(double mass_msun){
    double f = constants->kms/units->velocity; 
    return 50*pow(mass_msun/33. ,0.33) *f; /* Stone 1991 */
}
double runaway_velocity(double mass_code){
    const double runaway_frac_Ostar = 0.5; /* snl make runaway frac cfgparameters */
    const double runaway_frac_Bstar = 0.1;
    const double Ostar_mass = 20;
    const double Bstar_mass = 8;
    double mass_msun = mass_code * units->mass/constants->Msun;

    if( mass_msun > Ostar_mass &&
	cart_rand() < runaway_frac_Ostar){
	return  sample_exponential( runaway_mean_kick(mass_msun) );
    }
    if( mass_msun > Bstar_mass && mass_msun <= Ostar_mass &&
	cart_rand() < runaway_frac_Bstar ){
	return  sample_exponential( runaway_mean_kick(mass_msun) );
    }
    return 0;

}
void get_com_pos( int icell, double com_pos[], int level ){
    int i, idir, isign;
    int neighbors[num_neighbors];
    double dcom_pos[nDim];
    double sum[nDim];
    
    cell_center_position(icell,com_pos); /* start at cell position */

    cell_all_neighbors(icell,neighbors);
    for(idir=0;idir<nDim; idir++){ dcom_pos[idir]=0; sum[idir]=cell_gas_density(icell); } 
    for(i=0;i<num_neighbors; i++){
	idir = i/2;
	isign = 1 - (i % 2)*2; // even is +
	dcom_pos[idir] += isign*2*cell_size[level]*cell_gas_density(neighbors[i]);
	sum[idir] += cell_gas_density(neighbors[i]);
    }
    /* now assign COM position but stay within cell*/
    for(idir=0;idir<nDim; idir++){
	com_pos[idir] += (
	    dcom_pos[idir] > 0 ?
	    MIN( dcom_pos[idir]/sum[idir], cell_size[level]/2.01) :
	    MAX( dcom_pos[idir]/sum[idir],-cell_size[level]/2.01) );
    }
}

void grow_star_particle( int ipart, float dmass, int icell, int level) {
    int i;
    double add_mass;
    double new_density;
    double density_fraction, thermal_pressure;
    double sum_mass, pmass_orig;
    double com_pos[nDim];
	
    cart_assert( ipart < num_star_particles );
    add_mass = MIN( dmass, 0.667*cell_volume[cell_level(icell)]*cell_gas_density(icell) );
    
    /* add mass [at COM of neighbors] with cell momentum; */
    pmass_orig = particle_mass[ipart];
    sum_mass = pmass_orig + add_mass;
/* snl precision can be an issue */
/*     get_com_pos(icell, com_pos, level ); */
/*     cart_assert( cell_contains_position(icell,com_pos) );  */
/*     for ( i = 0; i < nDim; i++ ) { */
/*      particle_x[ipart][i] = (com_pos[i]*add_mass + particle_x[ipart][i]*pmass_orig)/sum_mass; */
/*     } */
    cart_assert( cell_contains_position(icell,particle_x[ipart]) ); 
    for ( i = 0; i < nDim; i++ ) {
	particle_v[ipart][i] = 
	    ( cell_momentum(icell,i) / cell_gas_density(icell) * add_mass +
	      particle_v[ipart][i]*pmass_orig ) /sum_mass;
    }
    particle_mass[ipart] += add_mass;
    star_initial_mass[ipart] += add_mass; /* this is used for feedback */
    
#ifdef ENRICHMENT
    if(sf_metallicity_floor>0.0 && cell_gas_metal_density_II(icell)<sf_metallicity_floor*constants->Zsun*cell_gas_density(icell)){
	cell_gas_metal_density_II(icell) =  sf_metallicity_floor*constants->Zsun*cell_gas_density(icell);
    }
    star_metallicity_II[ipart] = 
	( cell_gas_metal_density_II(icell) / cell_gas_density(icell) * add_mass +
	  star_metallicity_II[ipart] * pmass_orig) / sum_mass ;
#ifdef ENRICHMENT_SNIa
    star_metallicity_Ia[ipart] = 
	( cell_gas_metal_density_Ia(icell) / cell_gas_density(icell) * add_mass +
	  star_metallicity_Ia[ipart] * pmass_orig) / sum_mass;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
    
    /* adjust cell values */
    new_density = cell_gas_density(icell) - add_mass * cell_volume_inverse[level];
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
}
int create_star_particle( int icell, float mass, double pdt, int type ) {
	int i;
	int ipart;
	int id;
	int level;
	double uni[nDim];
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
#ifdef STAR_PARTICLE_TYPES
	if(starII_runaway_indicator == 1 && star_particle_type[ipart] == STAR_TYPE_STARII){
	    cart_rand_unit_vector(uni);
	    for ( i = 0; i < nDim; i++ ) {
		particle_v[ipart][i] += runaway_velocity(particle_mass[ipart])*uni[i];
	    }
	}
#endif

	particle_t[ipart] = tl[level];
	particle_dt[ipart] = pdt;

	star_tbirth[ipart] = tl[level];
	particle_mass[ipart] = true_mass;
	star_initial_mass[ipart] = true_mass;

#ifdef STAR_PARTICLE_TYPES
	if ( type == STAR_TYPE_NORMAL ) {
#endif /* STAR_PARTICLE_TYPES */

#ifdef ENRICHMENT
	star_metallicity_II[ipart] = cell_gas_metal_density_II(icell) / cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
	star_metallicity_Ia[ipart] = cell_gas_metal_density_Ia(icell) / cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */

#ifdef STAR_PARTICLE_TYPES
	} else {
#ifdef ENRICHMENT
    star_metallicity_II[ipart] = 0.0;
#ifdef ENRICHMENT_SNIa
	star_metallicity_Ia[ipart] = 0.0;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
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

	return ipart;
}

#endif /* HYDRO */
#endif /* PARTICLES && STAR_FORMATION */
