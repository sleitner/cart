#ifndef __FEEDBACK_KINETIC_H__
#define __FEEDBACK_KINETIC_H__
#ifdef HYDRO
void kfb_config_init();
void kfb_config_verify();
void distribute_momentum(double dp, int level, int icell, double dt);
void nonlocal_kicks(int level, int cell, double t_next, double dt );
#endif
#endif /* __FEEDBACK_STARII_H__ */
