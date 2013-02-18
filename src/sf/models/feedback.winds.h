#ifndef __FEEDBACK_WINDS_H__
#define __FEEDBACK_WINDS_H__

#ifdef STARFORM
void wind_config_init();
void wind_config_verify();
void wind_init();
void wind_setup(int level);
void stellar_wind_kick(int level,int cell,int ipart,double t_next);
#endif /* STARFORM */

#endif /* __FEEDBACK_WINDS_H__ */
