#include "config.h"

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "particle.h"
#include "starformation.h"
#include "units.h"
#include "onestarfits.h"
#include "form_star.starII.h"

#ifndef ENRICHMENT
#error "for onestarfits ENRICHMENT should be on"
#endif /* ENRICHMENT */

double agetau(double ini_mass_sol, double age_yr, double Zsol);
/* Raiteri et al. 1996 */
/* Mej=0.7682*pow(mass_sol,1.056) */
double OneStar_snII_Mejected_Fe(double mass_code){
    /* need units/cosmology set */
    return 2.802e-4 * pow(mass_code*units->mass/constants->Msun,1.864)/units->mass*constants->Msun;
}
double OneStar_snII_Mejected_Ox(double mass_code){
    /* need units/cosmology set */
    return 4.586e-4 * pow(mass_code*units->mass/constants->Msun,2.721)/units->mass*constants->Msun;
}

double OneStar_stellar_lifetime(double ini_mass_sol, double Zsol){  
    /* stellar lifetime in code units */
    /* need units/cosmology set */
    double logtstar, tstar;
    double a0,a1,a2;
    double logM, logZ;
    logM = log10(ini_mass_sol);
    /* for no enrichment use Z=-1.0, use Z=-3.0 floor*/
    /* logZ is in absolute metallicity*/
    logZ = log10( MAX( Zsol, 1e-2 )*constants->Zsun ); /* Z scales with oxygen abundance */ 
    a0 =  10.13 + 7.547e-2 *logZ - 8.084e-3 *(logZ*logZ);
    a1 = -4.424 - 7.939e-1 *logZ - 1.187e-1 *(logZ*logZ);
    a2 =  1.262 + 3.385e-1 *logZ + 5.417e-2 *(logZ*logZ);
    logtstar = a0 + a1*logM + a2*(logM*logM);
    tstar = pow(10,logtstar)*constants->yr/units->time;
    return tstar;
}

double agetau(double ini_mass_sol, double age_yr, double Zsol){
    double star_life;
    /* need units/cosmology set */
    star_life = OneStar_stellar_lifetime(ini_mass_sol, Zsol);
    return age_yr*constants->yr/units->time/star_life;
}
 

double OneStar_UV_fraction(double ini_mass_sol, double age_yr, double Zsol){
    const double sclZ=-2.0; 
    double logZsol, logZ_scl;
    double fmexp, fmlow, fmZ, ftZ, ft;
    double uvfraction, at;
#ifdef ENRICHMENT
    logZsol=log10(Zsol);
    logZsol=MIN( MAX(logZsol,-2), 0);  /* range of FSPS models -- extrapolation likely valid on either side.*/ 
    logZ_scl=logZsol/sclZ; /* Z=-2,0 -> 1,0 */
#else
    logZsol = -1;
    logZ_scl=0.5; 
#endif /* ENRICHMENT */
    if(ini_mass_sol > 5.0){
	    fmexp = 1 - exp( (5.0-ini_mass_sol)/27. );
	    fmlow = ini_mass_sol > 40 ? 1.0 : 50./(90-ini_mass_sol) ;
	    fmZ = 1 + logZ_scl * ( ini_mass_sol > 40 ? -.1 : 1)  ;
	    ftZ = 1 + 0.7*logZ_scl * agetau(ini_mass_sol, age_yr, Zsol);
	    ft = 0.55*(1.5 - agetau(ini_mass_sol, age_yr, Zsol));
	    ft = MAX( ft , 0.0 );
	    uvfraction =  fmlow * fmexp * fmZ * ftZ * ft;
	    uvfraction = MAX(0.0, MIN(1.0, uvfraction));
    }else{
	    uvfraction = 0;
    }
    return uvfraction;
}
double OneStar_ionizing_fraction(double ini_mass_sol, double age_yr, double Zsol){
    double fmlow, flolo, UV_frac, fion;
    fmlow = ini_mass_sol < 30 ? 0.25*((30-ini_mass_sol)/(30-8.)) : 0  ; 
    flolo = ini_mass_sol < 16 ? 1.3*((16-ini_mass_sol)/(16-8.)) : 0  ; 
    UV_frac = OneStar_UV_fraction(ini_mass_sol, age_yr, Zsol) ;
    fion = UV_frac * pow(10,-0.1-fmlow-flolo);
    fion = MAX(0.0, MIN(1.0, fion));
    return  fion;
}
double OneStar_Lbol_Lsun(double ini_mass_sol, double age_yr, double Zsol){
    double fmlow, fML, ftm, Lbol, tau; 
    tau = agetau(ini_mass_sol, age_yr, Zsol );
    if( tau > 3 || ini_mass_sol < 1.0 ){
	    Lbol = 0;
    }else{
	    fmlow = ini_mass_sol > 30 ? 1.0 : 1.0 - 0.07*pow( (30-(ini_mass_sol-8))/30. ,3) ;
	    fML = 3.9*pow( MAX(log10(ini_mass_sol),1.0) ,0.65);
	    ftm =1 - 0.1 * (ini_mass_sol-120)/120. * tau;
	    /* fit for > 8Msun stars only! */
	    Lbol = pow(10,fmlow*fML*ftm) ; /*in Lsun*/
    }
    cart_assert(Lbol>=0);
    return Lbol;
}

/* RT luminosity is in energy/restmassenergy/sec_codeunits */
const double Lsun_to_ergs = 3.826e33; /* Lsun ->ergs/s*/
#ifdef PARTICLES
#ifdef STARFORM
double OneStar_Lbol_RT(int ipart){
    double ini_mass_sol, age_yr, Zsol, restenergy, Lbol ; 
    /* need units/cosmology set */
    cart_error("careful -- RT wants an ionizing source function.");
    restenergy = particle_mass[ipart]*units->mass * constants->c*constants->c; /* erg */
    ini_mass_sol = star_initial_mass[ipart]*units->mass/constants->Msun;
    age_yr = (particle_t[ipart]-star_tbirth[ipart])*units->time/constants->yr ;
    Zsol = star_metallicity_II[ipart]/constants->Zsun;
    Lbol = OneStar_Lbol_Lsun(ini_mass_sol,age_yr,Zsol) * Lsun_to_ergs * units->time;
    return Lbol / restenergy ;
}
double OneStar_Lion_RT(int ipart){
    double ini_mass_sol, age_yr, Zsol, restenergy, Lion ; 
    /* need units/cosmology set */
    restenergy = particle_mass[ipart]*units->mass * constants->c*constants->c; /* erg */
    ini_mass_sol = star_initial_mass[ipart]*units->mass/constants->Msun;
    age_yr = (particle_t[ipart]-star_tbirth[ipart])*units->time/constants->yr ;
    Zsol = star_metallicity_II[ipart]/constants->Zsun;
    Lion = OneStar_Lbol_Lsun(ini_mass_sol,age_yr,Zsol) * Lsun_to_ergs * units->time
	    * OneStar_ionizing_fraction(ini_mass_sol,age_yr,Zsol);
    return Lion / restenergy ;
}
#endif /* STARFORM */
#endif /* PARTICLES */

double OneStar_Lbol_hydro(double ini_mass_sol, double age_yr, double Zsol){
    double hL =  OneStar_Lbol_Lsun(ini_mass_sol,age_yr,Zsol) * Lsun_to_ergs / units->energy * units->time;
    return hL;
}

double OneStar_wind_pdot_msunyrkms(double ini_mass_sol, double age_yr, double Zsol){ 
/* high-mass stellar wind momentum (\dot{M*}*vinf)*/
    const double zfac=0.8;
    const double yint=-6;
    const double m1=2.6;
    const double b1=4.35;
    const double b2=3.7; 
    const double m2=1.9;
    const double Lsplit=5.25;
    double ls,logZsol,logZ_scl,lmdotvw,mdotvw;
#ifdef ENRICHMENT
    logZsol=log10(Zsol);
    logZsol=MIN( MAX(logZsol,-2), 0);  /* range of FSPS models -- extrapolation likely valid on either side.*/
#else
    logZsol = -1;
#endif /* ENRICHMENT */
    logZ_scl=logZsol*zfac;
    
    ls = log10( OneStar_Lbol_Lsun(ini_mass_sol,age_yr,Zsol) );
    lmdotvw = ls < Lsplit ? m1*ls + (yint-b1*m1) : m2*ls + (yint-b2*m2);
    lmdotvw +=  logZ_scl;
    mdotvw = pow(10, lmdotvw);
    return mdotvw;
}
double OneStar_wind_pdot(double ini_mass_sol, double age_yr, double Zsol){ 
    double pdot =  OneStar_wind_pdot_msunyrkms(ini_mass_sol, age_yr, Zsol);
    /* need units/cosmology set */
    return pdot * 
	constants->Msun / constants->yr * constants->kms/
	(units->mass * units->velocity / units->time);
}
