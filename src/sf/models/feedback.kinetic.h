#ifndef __FEEDBACK_KINETIC_H__
#define __FEEDBACK_KINETIC_H__
#ifdef HYDRO
void kfb_config_init();
void kfb_config_verify();
void kfb_init();
void distribute_momentum(double dp, int level, int icell, double dt);
#endif
#endif /* __FEEDBACK_STARII_H__ */
