#ifndef __FEEDBACK_STARII_WINDS_H__
#define __FEEDBACK_STARII_WINDS_H__

#ifdef STARFORM
void starII_wind_config_init();
void starII_wind_config_verify();
void starII_stellar_wind_kick(int level, int icell, int ipart, double ini_mass_sol, double age_yr, double Zsol, double t_next);
#endif /* STARFORM */

#endif /* __FEEDBACK_STARII_WINDS_H__ */
