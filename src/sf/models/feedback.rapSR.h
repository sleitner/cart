#ifndef __FEEDBACK_RAPSR_H__
#define __FEEDBACK_RAPSR_H__

#ifdef STARFORM
#if defined(HYDRO) && defined(PARTICLES)
void rapSR_config_init();
void rapSR_config_verify();
void rapSR_setup(int level);

double rapSR_pdot(int ipart);
void rapSR_kick(int level, int icell, int ipart, double t_next);
#endif
#endif
#endif /* __FEEDBACK_RAPSR_H__ */
