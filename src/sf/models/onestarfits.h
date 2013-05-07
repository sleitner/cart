#ifndef __ONESTARFITS_H__
#define __ONESTARFITS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

double OneStar_snII_Mejected_Fe(double mass_code);
double OneStar_snII_Mejected_Ox(double mass_code);
double OneStar_stellar_lifetime(double ini_mass_sol, double Zsol);

double OneStar_wind_pdot(double ini_mass_sol, double age_yr, double Zsol);
double OneStar_wind_pdot_msunyrkms(double ini_mass_sol, double age_yr, double Zsol);

double OneStar_UV_fraction(double ini_mass_sol, double age_yr, double Zsol);
double OneStar_ionizing_fraction(double ini_mass_sol, double age_yr, double Zsol);

double OneStar_Lbol(double ini_mass_sol, double age_yr, double Zsol);
double OneStar_Lbol_Lsun(double ini_mass_sol, double age_yr, double Zsol);
double OneStar_Lbol_hydro(double ini_mass_sol, double age_yr, double Zsol);
#ifdef STARFORM 
#ifdef PARTICLES
double OneStar_Lbol_RT(int ipart);
double OneStar_Lion_RT(int ipart);
#endif
#endif
#ifdef DEBUG_SNL
void testonestar();
#endif

double agetau(double ini_mass_sol, double age_yr, double Zsol);

extern const double Lsun_to_ergs; /* Lsun ->ergs/s*/

#endif /* __ONESTARFITS_H__ */
