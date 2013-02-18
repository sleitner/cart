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

#ifdef HYDRO

void star_formation( int level, int time_multiplier );
void star_formation_normal( int level, int time_multiplier );

/* star formation parameters */
extern int sf_min_level;
extern float sf_metallicity_floor;
extern double sf_timescale;
extern double sf_min_stellar_particle_mass;

/* cluster formation parameters*/
void cluster_formation( int level, int time_multiplier);
int find_lowv_cluster(int icell);
extern int cluster_buildup_indicator;
extern double cluster_age_spread;     /* in yrs */
extern double cluster_min_expected_mass; /* in Msun*/
extern double fast_growth_probability;
extern double fast_growth_multiplier;

double cluster_age_spread_code; 
double cluster_min_expected_mass_code; 

#ifdef STAR_PARTICLE_TYPES
/* starII formation parameters*/
void starII_creation( double dmstarII, int icell, int level );
double msample_imf_highmass();
extern int starII_indicator; 
extern double starII_minimum_mass;
extern double starII_highmass_slope; 

extern double mfrac_starII, starII_avg_mass;
extern double sIIminpw, imfmaxpw;
double starII_avg_mass_code;
void setup_starII_formation(int level){
    starII_avg_mass_code = starII_avg_mass * constants->Msun / units->mass;
}

double msample_imf_highmass()
{
    /* Assume the high mass end is a power-law */
    /* should be IMF-dependent, but for now do Salpeter~Kroupa~Chabrier */
    return pow( cart_rand() * (imfmaxpw - sIIminpw)  + sIIminpw , 1/(starII_highmass_slope+1.0));
    /* invert CDF of P(m)=dN/dm */
    /* CDF = 1/Nhigh Int^M_Mmin m^-slope dm*/
    /* Nhigh= Int^Mmax_Mmin m^-slope dm*/
    /* CDF(m) = 1.35/1.35*(m^-1.35-Mmin^-1.35)/(Mmax^-1.35-Mmin^-1.35)*/
    /* m=CDF^-1(Puni)*/
}

void starII_creation( double dmstarII, int icell, int level ){
    int ipart;
    double mstargas_left, starII_mass, Pform_leftover ;
    /*  Form particles until alotted mass is used up. */

    mstargas_left = dmstarII ;
    while( mstargas_left > 0.0 ){

	starII_mass = msample_imf_highmass() *constants->Msun/units->mass;
	if( mstargas_left < starII_mass ){ /* last II mass forms stochastically. */

	    Pform_leftover = mstargas_left / starII_avg_mass_code; /* cannot be mass dependent: P(M*)=P(M*|IMF)*constant_wrt_M* */

	    if( cart_rand() < Pform_leftover ){ 
		ipart = create_star_particle( icell, starII_mass, dtl[level], STAR_TYPE_STARII ); 
	    }
	} else {

	    ipart = create_star_particle( icell, starII_mass, dtl[level], STAR_TYPE_STARII ); 

	}
	mstargas_left -= starII_mass;
    }
}
#endif /* STAR_PARTICLE_TYPES */

void setup_cluster_formation(int level){
    cluster_age_spread_code = cluster_age_spread * constants->yr/units->time ; 
    cluster_min_expected_mass_code = cluster_min_expected_mass * constants->Msun/units->mass ; 
}
int find_lowv_cluster(int icell){
    int ipart, ipart_store;
    double min_dv, dv, sage; 
    ipart = cell_particle_list[icell]; /* opt: could do oct instead of cell */
    min_dv=1e30;
    ipart_store=-1;
/* #ifdef DEBUG_SNL */
/*  if(ipart == NULL_PARTICLE){
	cart_debug("    no stars in cell %d",icell);
    }else{
 	cart_debug("    stars found in cell %d id%d",icell,ipart);*
    }
*/
/* #endif */
    while ( particle_is_star(ipart) 
#ifdef STAR_PARTICLE_TYPES
	    && (star_particle_type[ipart] == STAR_TYPE_NORMAL 
		|| star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH)
#endif
	){ /* if a cluster exists */
	/* if theres more than one seed then assign to star with lowest relative velocity; */
	sage = particle_t[ipart] - star_tbirth[ipart]; /* this can be negative (floats) */
	if (sage < cluster_age_spread_code){
	    dv =pow(particle_v[ipart][0]-cell_momentum(icell,0)/cell_gas_density(icell) , 2)+
		pow(particle_v[ipart][1]-cell_momentum(icell,1)/cell_gas_density(icell) , 2)+
		pow(particle_v[ipart][2]-cell_momentum(icell,2)/cell_gas_density(icell) , 2);
	    if( dv < min_dv ){ 
		min_dv = dv;
		ipart_store = ipart;
	    }
	}
/* #ifdef DEBUG_SNL  */
/* 	if(ipart_store == -1){  
 	    cart_debug(" found particles, but not young sage=%e spread=%e",sage,cluster_age_spread_code); 
 	}
*/
/*     cart_debug("min_dv %e",min_dv*units->velocity/constants->kms); */
/* #endif */
	ipart = particle_list_next[ipart];
    }
    return ipart_store;
}
 

/* creates cluster and builds until cluster age spread, samples starII; */
void cluster_formation( int level, int time_multiplier ) 
{
  int i, ipart;
  int icell;
  int num_level_cells;
  int *level_cells;
  double dmstar, dmcluster, dmstarII;
  float *sfr;

  if ( level < sf_min_level ) return;

  start_time( WORK_TIMER );
  select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
  sfr = cart_alloc(float,num_level_cells);
  star_formation_rate(level,num_level_cells,level_cells,sfr);

  
  for ( i = 0; i < num_level_cells; i++ ) { 
      if ( sfr[i] > 0 ){
	  icell = level_cells[i];
	  dmstar = sfr[i] * dtl[level]*time_multiplier * cell_volume[level];/* continous and deterministic */
	  dmcluster = dmstar;
	  dmstarII = 0;
#ifdef STAR_PARTICLE_TYPES
	  if( starII_indicator ) {
	      dmstarII = mfrac_starII*dmstar;
	      dmcluster = (1-mfrac_starII)*dmstar;
	      starII_creation( dmstarII, icell, level );
	  }
#endif /* STAR_PARTICLE_TYPES */
	  

	  ipart = find_lowv_cluster(icell); /* seeding or growing? */
	  if( ipart == -1 ){  /* no young clusters */ 
	      /* only seed clusters that you expect to grow to some mass*/
	      if ( sfr[i]*cluster_age_spread_code > cluster_min_expected_mass_code ){
#ifdef STAR_PARTICLE_TYPES
		  if(fast_growth_probability != 0 ){
		      if(cart_rand()<fast_growth_probability){
			  dmstarII *= fast_growth_multiplier;
			  dmcluster *= fast_growth_multiplier;
			  dmstar *= fast_growth_multiplier;
			  ipart = create_star_particle( icell, dmcluster, dtl[level], STAR_TYPE_FAST_GROWTH );
		      }
		  }else
#endif
		      {
			  ipart = create_star_particle( icell, dmcluster, dtl[level], STAR_TYPE_NORMAL ); 
		      }
		  star_initial_mass[ipart] += dmstarII; /* add back massII for feedback (mass loss is dealt with by starII) */
/* #ifdef DEBUG_SNL */
/*   double pos[nDim]; */
/* 		  cell_center_position(icell,pos); */
/* 		  cart_debug("seeding %d %d %e  %e %e %e    %f", */
/* 			     icell, ipart, */
/* 			     cell_gas_density(icell)*units->number_density/constants->cc, */
/* 			     pos[0],pos[1],pos[2], */
/* 			     dmstar*units->mass/constants->Msun */
/* 		      ); */
/* #endif */
		  
#ifdef BLASTWAVE_FEEDBACK
		  init_blastwave(icell);
#endif /* BLASTWAVE_FEEDBACK */
	      }
	  }else{
#ifdef STAR_PARTICLE_TYPES
	      if(star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH ){
		  dmstarII *= fast_growth_multiplier;
		  dmcluster *= fast_growth_multiplier;
		  dmstar *= fast_growth_multiplier;
	      }
#endif
	      grow_star_particle( ipart, dmcluster, icell, level);
	      star_initial_mass[ipart] +=  dmstarII; 
	  }
/* #ifdef DEBUG_SNL */
/* 	      cart_debug("growing id %d %e + %f",ipart, */
/* 			 particle_mass[ipart]*units->mass/constants->Msun, */
/* 			 dmstar*units->mass/constants->Msun); */
/* #endif */
	  
      }
  }
  
  cart_free( sfr );
  cart_free( level_cells );
  
  end_time( WORK_TIMER );
}




void star_formation_normal( int level, int time_multiplier )
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
		  mstar = MAX(mstar_min,mstar);
		else 
		  mstar = 0.0;
		  
#else  /* OLDSTYLE_SF_ALGORITHM */

		/* draw number of star formation events 0...\inf from poisson distribution */
		mstar = MAX( mstar_min, mstar );
		mstar *= (double)cart_rand_poisson( sfr[i]*cell_volume[level]*dt_eff/mstar );

#endif /* OLDSTYLE_SF_ALGORITHM */

		if ( mstar > 0.0 ) {

			/* create the new star */
			create_star_particle( icell, mstar, dtl[level], STAR_TYPE_NORMAL );

#ifdef LOG_STAR_CREATION      
			log_star_creation( icell, particle_mass[ipart], FILE_RECORD);
#endif

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

void star_formation( int level, int time_multiplier )
{

    if(cluster_buildup_indicator){
	cart_assert(time_multiplier == 1);
	setup_cluster_formation(level);
    }
#ifdef STAR_PARTICLE_TYPES
    if( starII_indicator ) {
	setup_starII_formation(level);
	cart_assert(cluster_buildup_indicator);
    } 
#endif /* STAR_PARTICLE_TYPES */
    if( cluster_buildup_indicator ){
	cluster_formation(level, time_multiplier);
    }else{
	star_formation_normal(level, time_multiplier);
    }
}

#endif /* STAR_FORMATION */
