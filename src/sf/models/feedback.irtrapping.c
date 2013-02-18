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

int Apply_AVK_tauIR = 0;
double tauIR_boost = 0;
double clump_dust_temp = 200;
double clump_survival_time = 5e6;
#ifdef STAR_PARTICLE_TYPES
extern double tdelay_popM_feedback;
#endif

void trapIR_config_init()
{
    control_parameter_add2(control_parameter_time,&clump_survival_time,"trapIR:clump-survival-time","clump_survival_time","survival time of clumps for AVK's tauIR.");

    control_parameter_add2(control_parameter_double,&clump_dust_temp,"trapIR:clump-dust-temp","clump_dust_temp","temperature of obscuring dust for tauIR (~constant above 125K).");

    control_parameter_add2(control_parameter_double,&tauIR_boost,"trapIR:boost","tauIR_boost","factor multiplying tauIR for RaP IR trapping.");

    control_parameter_add2(control_parameter_bool,&Apply_AVK_tauIR,"trapIR:AVK_model","Apply_AVK_tauIR","use Andrey's model for IR trapping based on matching clump and cluster mass functions.");
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
    const double eta2 = 1; /* clumping factor? */
    const double eps_clump = 0.2; /* star efficiency in clumps varies with \Sigma (see fall 2010). ; old 0.3 */
    const double mu_max = 1.0; /*  0.1 up to 1.0 ; old .1 */  
    const double Mclump_min = 100; 
    double Mclump_max;
    /* Sigma_cl - M_cl relation*/
    const double alpha0 = 0.4;
    const double Mflat1 = 3e4; /* for M>3e4 (Fig1) */
    const double alpha1 = 0.0; 
    double alpha ;
    const double beta = 1.7; /* slope of the clump mass function */
    double ab,b2;
    double C_R;
    double tauIR;

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

void masslum_from_star0(int level, int icell, double *Msol, double *LUV_ergis){
    int ipart;
    double mstar, lstar; 

    mstar=0;
#ifndef RADIATIVE_TRANSFER
    lstar=0;
#endif
    ipart = cell_particle_list[icell];
    while ( ipart != NULL_PARTICLE ) {
 	if ( particle_is_star(ipart) ){
	    if( (particle_t[ipart]-star_tbirth[ipart])*units->time 
		< constants->yr*clump_survival_time ){
		mstar += particle_mass[ipart];
#ifndef RADIATIVE_TRANSFER
		lstar += sf_feedback->rt_source(ipart)*particle_mass[ipart]; 
#endif
	    }
	}
	ipart = particle_list_next[ipart];
    }
#ifdef RADIATIVE_TRANSFER
    /* uses cic (rtAfterAssignDensity adds inv_volume) */
    (*LUV_ergis) = cell_rt_source(icell)*cell_volume[level]  * factor_rt_to_ergis ; 
#else
    (*LUV_ergis) = lstar * factor_rt_to_ergis; /* uses ngp */             
#endif
    (*Msol) =  mstar*units->mass/constants->Msun;
}


void cell_trapIR(int level, int icell, double t_next, double dt){
    double Zsol;
    double Lbol_ergis, LUV_ergis;
    double Mcell_sun;
    double tauIR;
    double dp;
    
    if( tauIR_boost > 0 ){

#ifdef ENRICHMENT
	Zsol = cell_gas_metal_density(icell)/(cell_gas_density(icell)*constants->Zsun);
#else
	Zsol = 1.0;
#endif
	masslum_from_star0(level, icell, &Mcell_sun, &LUV_ergis);
	Lbol_ergis = LUV_ergis*6; /* similar to rapSR_luminosity_code*/
	if(LUV_ergis>0){
/* #ifdef DEBUG_SNL */
/* 	    cart_debug("mass,lum from star0 %e %e",Mcell_sun, LUV_ergis);  */
/* #endif */
	    
	    tauIR = tauIR_boost *  Kappa_IR(Zsol) * constants->XH
#ifdef RADIATIVE_TRANSFER
		* (cell_HI_density(icell)+2.0*cell_H2_density(icell))
#else
		* cell_gas_density(icell)
#endif /* RADIATIVE_TRANSFER */
		* Zsol * cell_sobolev_length(icell);
	    
	    if(Apply_AVK_tauIR)
		{
		    tauIR += tauIR_boost * AVK_tauIR(Mcell_sun, Zsol) ;
/* #ifdef DEBUG_SNL */
/* 		double tauIR2; */
/* 		tauIR2 = tauIR_boost *  Kappa_IR(Zsol) * constants->XH  */
/* #ifdef RADIATIVE_TRANSFER */
/* 		  * (cell_HI_density(icell)+2.0*cell_H2_density(icell))  */
/* #else */
/* 		  * cell_gas_density(icell) */
/* #endif */
/* 		  * Zsol * cell_sobolev_length(icell); */
/* 		cart_debug("AVK trapped tauIR=%e cell tauIR=%e",tauIR, tauIR2); */
/* #endif /\* DEBUG_SNL *\/ */
		}
	    dp = LUV_ergis * tauIR / constants->c * (dt*units->time) 
		/ (units->velocity * units->mass);
/*  	    cart_debug("LIR %e dpIR %e %e",LUV_ergis, dp/tauIR, dp);  */
	    distribute_momentum(dp, level, icell, dt);
	}
    }
}

#endif /* STAR_FORMATION */
