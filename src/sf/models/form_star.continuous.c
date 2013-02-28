#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include <math.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "particle.h"
extern int continuous_starformation_indicator;   
#include "starformation.h"
#include "starformation_recipe.h"
#include "rand.h"
#include "rt.h"
#include "tree.h"
#include "units.h"

#include "form_star.continuous.h"
#include "form_star.starII.h"

int find_lowv_cluster(int icell);
double fast_growth_probability=0;            /* how often does efficient star formation happen? */
double fast_growth_multiplier=10;            /* how quickly does the efficient star formation happen? */
double cluster_age_spread=15.0e6;            /* in yrs; age above which "cluster" formation ceases */
double cluster_min_expected_mass=100;        /* do not form a cluster particle if it is expected to grow to less than this mass */

double cluster_age_spread_code; 
double cluster_min_expected_mass_code; 
int starII_indicator;

void continuous_config_init()
{
//    if(!(continuous_starformation_indicator)) return;
    control_parameter_add3(control_parameter_time, &cluster_age_spread,"cluster:age-spread","cluster_age_spread","cluster_age_spread","timescale over which star particles representing clusters are allowed to grow.");
    control_parameter_add2(control_parameter_double, &cluster_min_expected_mass,"cluster:min-expected-mass","cluster_min_expected_mass","the minimum mass expected from sfr*cluster_age_spread allowed to seed a cluster.");
    
    control_parameter_add2(control_parameter_double, &fast_growth_multiplier, "cluster:fast-growth-multiplier", "fast_growth_multiplier", "multiplier for fast growth particle star formation rates");
    control_parameter_add2(control_parameter_double, &fast_growth_probability, "cluster:fast-growth-probability", "fast_growth_probability", "probability of a fast growth particle");

    control_parameter_add2(control_parameter_bool, &starII_indicator, "starII:indicator", "starII_indicator", "turn on starII star-formation");
    starII_config_init();
}

extern double sf_sampling_timescale ;        /* in yrs; used to be called dtmin_SF, also in HART */
void continuous_config_verify()
{
//   if(!(continuous_starformation_indicator)) return;

    VERIFY(cluster:fast-growth-multiplier, fast_growth_multiplier >= 0);//cluster:fast-growth-multiplier
    VERIFY(cluster:fast-growth-probability, fast_growth_probability >= 0 && fast_growth_probability <=1 );
    
    VERIFY(sf:sampling-timescale, sf_sampling_timescale == 0 );

    VERIFY(cluster:age-spread, cluster_age_spread >1.0e6 );
    VERIFY(cluster:min-expected-mass, cluster_min_expected_mass > 0.0 );

    VERIFY(starII:indicator,starII_indicator==1 || starII_indicator==0);
    starII_config_verify();
}
void continuous_init(){
    if(!(continuous_starformation_indicator)) return;
    starII_init();
}
void continuous_setup(int level){
    if(!(continuous_starformation_indicator)) return;
    cluster_age_spread_code = cluster_age_spread * constants->yr/units->time ; 
    cluster_min_expected_mass_code = cluster_min_expected_mass * constants->Msun/units->mass ; 
    starII_setup(level);
}
int find_lowv_cluster(int icell){
    int ipart, ipart_store;
    double min_dv, dv, sage; 
    ipart = cell_particle_list[icell]; /* opt: could do oct instead of cell */
    min_dv=1e30;
    ipart_store=-1;
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
	ipart = particle_list_next[ipart];
    }
    return ipart_store;
}
 
extern double mfrac_starII ; /* this is given by IMF */
void continuous_star_formation( int level, int icell, double dtl, double dt_eff, float sfr) {
/* creates cluster and builds until cluster age spread, samples starII; */
    int i, ipart;
    int num_level_cells;
    int *level_cells;
    float dmstar;
    int star_type = STAR_TYPE_NORMAL;
    
    if(!(continuous_starformation_indicator)) return;
    if ( sfr <= 0 ) return;
    
    dmstar = sfr * dt_eff * cell_volume[level];
    
    ipart = find_lowv_cluster(icell); /* seeding or growing? */
    if( ipart == -1 ){  
        /* SEED -- there are no young clusters  -----------------*/ 
        
        /* seed clusters that you expect to grow to some mass -- amounts to a density threshold */
        if ( sfr*cluster_age_spread_code < cluster_min_expected_mass_code )  return; 
	
#ifdef STAR_PARTICLE_TYPES
        if(fast_growth_probability != 0 ){
            if(cart_rand()<fast_growth_probability){ 
                dmstar *= fast_growth_multiplier; 
                star_type = STAR_TYPE_FAST_GROWTH;
            }
        }
        if( starII_indicator ) { starII_creation( mfrac_starII*dmstar, icell, level, dtl );}
#endif /* STAR_PARTICLE_TYPES */
        ipart = create_star_particle( icell, (1-mfrac_starII)*dmstar, dtl, star_type ); 
            
    }else{ 
        /* GROW -- there is a young cluster present ------------*/
       
#ifdef STAR_PARTICLE_TYPES
        if(star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH ){ dmstar *= fast_growth_multiplier;}
        if( starII_indicator ) { starII_creation( mfrac_starII*dmstar, icell, level, dtl );}
#endif
        grow_star_particle( ipart, (1-mfrac_starII)*dmstar, icell, level);
    }
    if(ipart != -1)
        star_initial_mass[ipart] += mfrac_starII*dmstar; /* add back massII for feedback (mass loss is dealt with by starII) */
  
}

#endif /* defined(HYDRO) && defined(STAR_FORMATION) */

