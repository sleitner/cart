#include "config.h"
#ifdef STAR_FORMATION

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "imf.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "tree.h"
#include "units.h"
#include "rand.h"

#include "form_star.runaway-starII.h"
#include "onestarfits.h"

/* STARII related */
double starII_highmass_slope=-2.35;          /* IMF slope used for sampling starII masses */
double starII_minimum_mass=8.0;              /* in Msun; mass above which individual stars are sampled */
int starII_runaway_indicator=0;              /* are we letting individual massive stars runaway? */

double starII_avg_mass_code;

double mfrac_starII=0, starII_avg_mass, m_imf_tot;
double sIIminpw, imfmaxpw;

double tdelay_popM_feedback;
void starII_config_init()
{
    control_parameter_add2(control_parameter_double,&starII_highmass_slope,"starII:highmass-slope","starII_highmass_slope","IMF slope used for sampling starII masses.");
    control_parameter_add2(control_parameter_double,&starII_minimum_mass,"starII:minimum-mass","starII_minimum_mass","the minimum mass of 'virtual' starII particles in Msun.");
    
    control_parameter_add2(control_parameter_bool, &starII_runaway_indicator, "starII:runaway-indicator", "starII_runaway_indicator", "turn on runaway starIIs");
    starII_runaway_config_init();
}

#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)
void check_fbdefs_compatible()
{
#ifdef SF_FEEDBACK
    const char *feedback_external_name = to_string(SF_FEEDBACK);
#else
    const char *feedback_external_name = "";
#endif
    if(strcmp("<popM-starII>",feedback_external_name)!=0){
        cart_error("SF_FORMSTAR includes STARII formation then SF_FEEDBACK must be a -starII variant : %s", feedback_external_name);
    }
}
void starII_config_verify()
{
#ifndef STAR_PARTICLE_TYPES
    cart_error("STAR_PARTICLE_TYPES must be defined for starII formation");
#endif /* STAR_PARTICLE_TYPES */
	
    VERIFY(starII:runaway-indicator,starII_runaway_indicator==1 || starII_runaway_indicator==0);
    VERIFY(starII:highmass-slope,starII_highmass_slope);
    VERIFY(starII:minimum-mass, starII_minimum_mass >1.0 );
    if(starII_runaway_indicator == 1){
	check_fbdefs_compatible();
    }
    starII_runaway_config_verify();
}   

void starII_init()
{
    double lowmass_avg_mass;
    double Nfrac_lowmass, Nfrac_starII;
    double N_imf_tot, N_imf_lowmass, N_imf_starII ;
    double m_imf_lowmass, m_imf_starII; 
    double mfrac_lowmass;
    
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

void starII_setup(int level){
    starII_avg_mass_code = starII_avg_mass * constants->Msun / units->mass;
    tdelay_popM_feedback = OneStar_stellar_lifetime(starII_minimum_mass, 1.0); /* need units */
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

void starII_creation( double dmstarII, int icell, int level, double dtl ){
    int ipart, i;
    double mstargas_left, starII_mass, Pform_leftover ;
    double vadd[nDim];
    /*  Form particles until alotted mass is used up. */
    mstargas_left = dmstarII ;
    while( mstargas_left > 0.0 ){
	starII_mass = msample_imf_highmass() *constants->Msun/units->mass;
	if( mstargas_left < starII_mass ){ /* last II mass forms stochastically. */
	    Pform_leftover = mstargas_left / starII_avg_mass_code; /* cannot be mass dependent: P(M*)=P(M*|IMF)*constant_wrt_M* */
	    if( cart_rand() < Pform_leftover ){ 
		ipart = create_star_particle( icell, starII_mass, dtl, STAR_TYPE_STARII ); 
	    }
	} else {
	    ipart = create_star_particle( icell, starII_mass, dtl, STAR_TYPE_STARII ); 
	}
        
        if(starII_runaway_indicator){
	    if(starII_runaway_velocity(starII_mass, vadd)){
		for ( i = 0; i < nDim; i++ ) {
		    particle_v[ipart][i] += vadd[i];
		}
	    }
	}
        mstargas_left -= starII_mass;
    }
}
#endif /* STAR_FORMATION */

