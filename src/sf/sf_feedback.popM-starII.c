#include "config.h"
#ifdef STAR_PARTICLE_TYPES
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "tree.h"
#include "units.h"
#include "starformation_feedback.h"
#include "models/onestarfits.h"

#include "models/feedback.snII.h"
#include "models/feedback.snIa.h"
#include "models/feedback.ml.h"
#include "models/feedback.winds.h"
#include "models/feedback.rapSR.h"
#include "models/feedback.rad.h"

#include "models/feedback.starII-snII.h"
#include "models/feedback.starII-winds.h"
#include "models/feedback.starII-rapSR.h"
#include "models/feedback.starII-rad.h"

#include "models/feedback.irtrapping.h"
#include "models/feedback.kinetic.h"

#include "times.h"
#include "timing.h"
#include "starformation.h"
#include "particle.h"

extern double starII_minimum_mass;
extern double tdelay_popM_feedback;

extern double rapSR_boost;
extern double wind_momentum_boost;
extern double starII_rapSR_boost;
extern double starII_wind_momentum_boost;

extern double tauIR_boost;

void sfb_config_init()
{
  snII_config_init();
  snIa_config_init();
  ml_snl2012_config_init();
  wind_config_init();
  rapSR_config_init();

  starII_wind_config_init();
  starII_rapSR_config_init();
  starII_explosion_config_init();

  trapIR_config_init();
  kfb_config_init();
}



#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)
void check_fsdefs_compatible()
{
#ifdef SF_FORMSTAR
    const char *formstar_external_name = to_string(SF_FORMSTAR);
#else
    const char *formstar_external_name = "";
#endif
    if(strcmp("<continuous>",formstar_external_name)!=0){
        cart_error("SF_FORMSTAR needs to be <continous> for SF_FEEDBACK -starII variants");
    }
}
void sfb_config_verify()
{
  check_fsdefs_compatible();

  snII_config_verify();
  snIa_config_verify();
  ml_snl2012_config_verify();
  wind_config_verify();
  rapSR_config_verify();

  starII_wind_config_verify();
  starII_rapSR_config_verify();
  starII_explosion_config_verify();

  trapIR_config_verify();
  kfb_config_verify();
}


void sfb_init()
{
  snII_init();
  snIa_init();
  ml_init();
  wind_init();
  rapSR_init(); 

/*   starII_wind_init(); */
/*   starII_rapSR_init(); */
  starII_explosion_feedback_init();

/*   trapIR_init();  */
}

void sfb_setup(int level)
{
    snII_setup(level);
    snIa_setup(level);
    rad_setup(level);
    ml_setup(level);
    wind_setup(level);
    rapSR_setup(level);

/*     starII_wind_setup(level); */
    starII_rapSR_setup(level);
    starII_explosion_setup(level);

    trapIR_setup(level);
/*     kfb_setup(level); */
}


#if defined(HYDRO) && defined(PARTICLES)
void sfb_hydro_feedback(int level, int cell, int ipart, double t_next )
{
    double ini_mass_sol, Zsol ;
    float star_age = particle_t[ipart] - star_tbirth[ipart] ; /* this can be negative (float) */
    double age_yr = star_age * units->time/constants->yr;
    Zsol = star_metallicity_II[ipart]/constants->Zsun;

    if( star_particle_type[ipart] ==  STAR_TYPE_NORMAL || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH ){ 

 	if( star_age > tdelay_popM_feedback ){  /* turn on popM fb after starII feedback is done */
	    /* not exactly right -- big dt misses the start when t_next>lifetime */
	    snII_thermal_feedback(level,cell,ipart,t_next);
	    ml_feedback(level,cell,ipart,t_next);
	    
 	    snII_kinetic_feedback(level,cell,ipart,t_next); 
	    if(rapSR_boost > 0){
		rapSR_kick(level,cell,ipart,t_next);    
	    }
	    if(wind_momentum_boost > 0){
		stellar_wind_kick(level,cell,ipart,t_next);
	    }
	}
	snIa_thermal_feedback(level,cell,ipart,t_next); 

    } else if ( star_particle_type[ipart] == STAR_TYPE_STARII ){

	ini_mass_sol = star_initial_mass[ipart]*units->mass/constants->Msun;

	if( star_age > OneStar_stellar_lifetime(ini_mass_sol, Zsol) ){ 

	    starII_explosion_mass(level, cell, ipart);
 	    starII_explosion_kicks(level, cell, ipart); 
	    starII_explosion_thermal(level, cell, ipart);  
	}else{
	    
	    if(cell_gas_density(cell) != cell_gas_density(cell)){
		cart_error("bad gas density %d %e",cell,cell_gas_density(cell) );
	    }

	    if(starII_rapSR_boost > 0){
		starII_rapSR_kick(level, cell, ipart,ini_mass_sol,age_yr,Zsol, t_next);    
	    }
	    if(starII_wind_momentum_boost > 0){
		starII_stellar_wind_kick(level, cell, ipart,ini_mass_sol,age_yr,Zsol, t_next); 
	    }
	    /* RaP longrange is cell-by-cell */
	} 
    }

}

int sfb_destroy_star_particle(int level,int icell,int ipart)
{
    double star_age, ini_mass_sol, Zsol;
    star_age = particle_t[ipart] - star_tbirth[ipart] ;
    ini_mass_sol = star_initial_mass[ipart]*units->mass/constants->Msun;
    Zsol = star_metallicity_II[ipart]/constants->Zsun;
    if ( star_particle_type[ipart] == STAR_TYPE_STARII &&
	 star_age > OneStar_stellar_lifetime(ini_mass_sol, Zsol)
	){
	return -1;
    }else{
	return 1;
    }
    
}


extern double sf_min_gas_number_density;
void sfb_hydro_feedback_cell(int level, int cell, double t_next, double dt )
{
    if( cell_gas_density(cell)*units->number_density*constants->XH 
	> sf_min_gas_number_density){ /* note this parameter is active for any SFP */
	if(tauIR_boost>0){
	    cell_trapIR(level, cell, t_next, dt); 
	} 
    }
}
#endif /* HYDRO && PARTICLES */

struct StellarFeedbackParticle sf_feedback_particle_internal =
{
    "popM-starII",
    sfb_hydro_feedback,
    rad_luminosity_popM_ionizingstarII0,
    NULL,
    sfb_config_init,
    sfb_config_verify,
    sfb_init,
    sfb_setup,
    sfb_destroy_star_particle
};

struct StellarFeedbackCell sf_feedback_cell_internal = 
{
    sfb_hydro_feedback_cell, 
    nonlocal_kicks
};

#endif /* STARFORM */
#endif /* STAR_PARTICLE_TYPES */

