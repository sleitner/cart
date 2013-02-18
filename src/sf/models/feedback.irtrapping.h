#ifndef __FEEDBACK_IRTRAPPING_H__
#define __FEEDBACK_IRTRAPPING_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void trapIR_config_init();
void trapIR_config_verify();
void trapIR_setup(int level);
void cell_trapIR(int level, int icell, double t_next, double dt);
double AVK_tauIR(double Mcell_sun, double Zsol) ;


#endif /* __FEEDBACK_IRTRAPPING_H__ */ 
