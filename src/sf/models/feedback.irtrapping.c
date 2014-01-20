#include "config.h"
#ifdef STAR_FORMATION
#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "imf.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "tree.h"
#include "units.h"

#include "onestarfits.h"
#include "feedback.kinetic.h"
#include "feedback.starII-rapSR.h"
#include "feedback.rapSR.h"

int Apply_AVK_tauIR = 0;
double tauIR_boost = 0;
double clump_dust_temp = 200;
double clump_survival_time = 5e6;

void trapIR_config_init()
{
    control_parameter_add2(control_parameter_time,&clump_survival_time,"trapIR:clump-survival-time","clump_survival_time","survival time of clumps for AVK's tauIR.");

    control_parameter_add2(control_parameter_double,&clump_dust_temp,"trapIR:clump-dust-temp","clump_dust_temp","temperature of obscuring dust for tauIR (~constant above 125K).");

    control_parameter_add2(control_parameter_double,&tauIR_boost,"trapIR:boost","tauIR_boost","factor multiplying tauIR for RaP IR trapping.");

    control_parameter_add2(control_parameter_bool,&Apply_AVK_tauIR,"trapIR:AVK_model","Apply_AVK_tauIR","use Andrey's model for IR trapping based on matching clump and cluster mass functions, if density is over the SF density threshold.");
}

void trapIR_config_verify()
{
    VERIFY(trapIR:clump-dust-temp, clump_dust_temp>0);
    VERIFY(trapIR:clump-survival-time, clump_survival_time>1000 && clump_survival_time<5e7);
    VERIFY(trapIR:boost, tauIR_boost >=0 );
}

static double factor_Kappa_IR;
static double factor_rt_to_ergis ;
void trapIR_setup(int level)
{
	if(clump_dust_temp < 125){
	/* Low temp  (<125K) Rosseland mean opacity from Semenov 2003 (Krumholz & Thompson 2012) */
	/* snl: added + 0.1 to be closer to median of lines */
		factor_Kappa_IR = tauIR_boost  
			* (pow(10,-1.5)*clump_dust_temp*clump_dust_temp/100. + 0.1)
			* units->mass/(units->length*units->length) ;
		
    }else{
	    /* High temp (>100K) Rosseland mean opacity from Semenov 2003 (Hopkins 2011) */
	    factor_Kappa_IR = tauIR_boost * 5 *units->mass/(units->length*units->length)  ; //5cm^2/gram  */
    }
    factor_rt_to_ergis =  units->mass * constants->c*constants->c / units->time ; 
}

double Kappa_IR(double Zsol){
    return factor_Kappa_IR*Zsol;
}

double AVK_tauIR(double Msol, double Zsol ){
    /* Andrey's scheme for distributing trapped IR radiation to clumps of different Sigma*/
	const float eta2 = 1; /* clumping factor? */
    const float eps_clump = 0.2; /* star efficiency in clumps varies with \Sigma (see fall 2010). ; old 0.3 */
    const float mu_max = 1.0; /*  0.1 up to 1.0 ; old .1 */  
    const float Mclump_min = 100; 
    float Mclump_max;
    /* Sigma_cl - M_cl relation*/
    const float alpha0 = 0.4;
    const float Mflat1 = 3e4; /* for M>3e4 (Fig1) */
    const float alpha1 = 0.0; 
    float alpha ;
    const float beta = 1.7; /* slope of the clump mass function */
    float ab,b2;
    float C_R;
    float tauIR;

     alpha = alpha1; 
     Mclump_max = mu_max * Msol / eps_clump ; 

/*     Rh = C_R * pow(Mclump, alpha ); */
     C_R = 2.5 * constants->pc * pow( Mflat1*constants->Msun, -alpha) ; 
/*     Mclump = Msol/eps_clump; */
/*     Sigma_cl = (1-eps_clump)*Mclump/(2*M_PI*Rh*Rh); */
     
     ab = 3-2*alpha-beta;
     b2 = 2-beta;
    
     tauIR = eta2 * Kappa_IR(Zsol) *units->length*units->length/units->mass 
	 * (1-eps_clump)/(2*M_PI*C_R*C_R)
	 * b2/ab * pow(mu_max/eps_clump,1-2*alpha)
	 * (1 - pow(Mclump_min/Mclump_max, ab)) / (1 - pow(Mclump_min/Mclump_max, b2) )
	 * pow(Msol * constants->Msun, 1-2*alpha) ;  

     tauIR = MIN(200.0, tauIR);
     return tauIR;
}

void masspdot_from_cell(int level, int icell, float *Msol_cell, float *pdot_cell){
    int ipart;
    double mstar_sum, pdot_sum; 
    double ini_mass_sol, Zsol, age_yr; 
    mstar_sum=0;
    pdot_sum=0;
    ipart = cell_particle_list[icell];
    while ( ipart != NULL_PARTICLE ) {
	    if ( particle_is_star(ipart) ){
		    age_yr = (particle_t[ipart]-star_tbirth[ipart])*units->time / constants->yr;
		    if( age_yr < clump_survival_time ){
			    mstar_sum += particle_mass[ipart];
#ifdef STAR_PARTICLE_TYPES
			    if(star_particle_type[ipart] == STAR_TYPE_NORMAL 
			       || star_particle_type[ipart] == STAR_TYPE_FAST_GROWTH ){
				    pdot_sum += rapSR_pdot(ipart);
			    }else if(star_particle_type[ipart] == STAR_TYPE_STARII){
				    ini_mass_sol = star_initial_mass[ipart]*units->mass/constants->Msun;
				    Zsol = star_metallicity_II[ipart]/constants->Zsun;
				    pdot_sum += starII_rapSR_pdot( ini_mass_sol, age_yr, Zsol) ;
			    }else{
				    cart_error("star type  needs an associated pdot");
			    }
#else
			    pdot_sum += rapSR_pdot(ipart);
#endif
		    }
	    }
	    ipart = particle_list_next[ipart];
    }
    (*pdot_cell) = pdot_sum; /* uses ngp */             
    (*Msol_cell) = mstar_sum*units->mass/constants->Msun;
}
/* #ifndef RADIATIVE_TRANSFER */
/* 		pdot += sf_feedback_particle->rt_source(ipart)*particle_mass[ipart]* factor_rt_to_ergis;  */
/* #endif */
/* #ifdef RADIATIVE_TRANSFER */
/*     /\* uses cic (rtAfterAssignDensity adds inv_volume) *\/ */
/*     (*pdot_cell) = cell_rt_source(icell)*cell_volume[level]  */
/* 	* 6 /\* Lbol ~ 6*Lion *\/ */
/* 	* factor_rt_to_ergis / constants->c  */
/* 	* (units->time)/ (units->velocity * units->mass); */
/* #endif */
extern double sf_min_gas_number_density;
void cell_trapIR(int level, int icell, double t_next, double dt){
    double Zsol;
    float pdot_cell;
    float Msol_cell;
    double tauIR;
    double dp;
#ifndef ENRICHMENT
    cart_error("ENRICHMENT is required for tauIR model to function");
#endif
    if( tauIR_boost > 0 &&
        cell_gas_density(icell)*units->number_density*constants->XH > sf_min_gas_number_density
        ){
#ifdef ENRICHMENT
	    Zsol = cell_gas_metal_density(icell)/(cell_gas_density(icell)*constants->Zsun);
#else
	    Zsol = 0;
#endif /* ENRICHMENT */
	    masspdot_from_cell(level, icell, &Msol_cell, &pdot_cell);
	    //Lbol_ergis = LUV_ergis*6; /* similar to rapSR_luminosity_code*/
	    if(pdot_cell>0){
		    tauIR = tauIR_boost *  Kappa_IR(Zsol) 
#ifdef RADIATIVE_TRANSFER
			    * (cell_HI_density(icell)+2.0*cell_H2_density(icell))/constants->XH
#else
			    * cell_gas_density(icell)
#endif /* RADIATIVE_TRANSFER */
			    * Zsol * cell_sobolev_length(icell);
		    if(Apply_AVK_tauIR){
			    tauIR += tauIR_boost * AVK_tauIR(Msol_cell, Zsol) ;
		    }
		    dp = pdot_cell * tauIR * dt ;
		    distribute_momentum(dp, level, icell, dt);
	    }
    }
}

#endif /* STAR_FORMATION */
