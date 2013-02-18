#ifndef __FEEDBACK_STARII_RAPSR_H__
#define __FEEDBACK_STARII_RAPSR_H__
#ifdef STARFORM
#if defined(HYDRO) && defined(PARTICLES)
void starII_rapSR_config_init();
void starII_rapSR_config_verify();
void starII_rapSR_setup();

void starII_rapSR_kick(int level, int icell, int ipart, double ini_mass_sol, double age_yr, double Zsol, double t_next);
#endif
#endif
#endif /* __FEEDBACK_STARII_RAPSR_H__ */
